#-Begin-preamble-------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 7.5.0
#
#
#     Main author:          Giovanni IADAROLA
#                           BE-ABP Group
#                           CERN
#                           CH-1211 GENEVA 23
#                           SWITZERLAND
#                           giovanni.iadarola@cern.ch
#
#     Contributors:         Eleonora Belli
#                           Philipp Dijkstal
#                           Lotta Mether
#                           Annalisa Romano
#                           Giovanni Rumolo
#
#
#     Copyright  CERN,  Geneva  2011  -  Copyright  and  any   other
#     appropriate  legal  protection  of  this  computer program and
#     associated documentation reserved  in  all  countries  of  the
#     world.
#
#     Organizations collaborating with CERN may receive this program
#     and documentation freely and without charge.
#
#     CERN undertakes no obligation  for  the  maintenance  of  this
#     program,  nor responsibility for its correctness,  and accepts
#     no liability whatsoever resulting from its use.
#
#     Program  and documentation are provided solely for the use  of
#     the organization to which they are distributed.
#
#     This program  may  not  be  copied  or  otherwise  distributed
#     without  permission. This message must be retained on this and
#     any other authorized copies.
#
#     The material cannot be sold. CERN should be  given  credit  in
#     all references.
#
#-End-preamble---------------------------------------------------------

import numpy as np
from numpy import sqrt, exp
from numpy.random import rand
import electron_emission as ee

def yield_fun2(E, costheta, Emax, del_max, R0, E0, s, flag_costheta_delta_scale=True, flag_costheta_Emax_shift=True):

    if flag_costheta_delta_scale:
        del_max_tilde=del_max*exp(0.5*(1.-costheta))
    else:
        del_max_tilde=del_max

    if flag_costheta_Emax_shift:
        E_max_tilde=Emax*(1.+0.7*(1.-costheta))
    else:
        E_max_tilde=Emax

    x=E/E_max_tilde

    true_sec=del_max_tilde*(s*x)/(s-1.+x**s)
    reflected=R0*((sqrt(E)-sqrt(E+E0))/(sqrt(E)+sqrt(E+E0)))**2.

    delta=true_sec+reflected
    ref_frac=0.*delta
    mask_non_zero=(delta>0)
    ref_frac[mask_non_zero]=reflected[mask_non_zero]/delta[mask_non_zero]

    return delta, ref_frac


class SEY_model_ECLOUD:
    def __init__(   
                    self, Emax,del_max,R0,
                    E_th=None, sigmafit=None, mufit=None, 
                    switch_no_increase_energy=0, thresh_low_energy=None,secondary_angle_distribution=None, 
                    E0=150., s=1.35, flag_costheta_delta_scale=True, flag_costheta_Emax_shift=True
                ):

        self.E_th = E_th
        self.sigmafit = sigmafit
        self.mufit = mufit
        self.switch_no_increase_energy = switch_no_increase_energy
        self.thresh_low_energy = thresh_low_energy
        self.secondary_angle_distribution = secondary_angle_distribution

        if secondary_angle_distribution is not None:
            import electron_emission
            self.angle_dist_func = electron_emission.get_angle_dist_func(secondary_angle_distribution)
        else:
            self.angle_dist_func = None

        self.Emax = Emax
        self.del_max = del_max
        self.R0 = R0
        self.E0 = E0
        self.s = s
        self.flag_costheta_delta_scale = flag_costheta_delta_scale
        self.flag_costheta_Emax_shift = flag_costheta_Emax_shift

        print 'Secondary emission model: ECLOUD E0=%.4f s=%.4f' % (self.E0, self.s)

    def SEY_process(self,nel_impact,E_impact_eV, costheta_impact, i_impact):
            
        yiel, ref_frac = yield_fun2(
            E_impact_eV,costheta_impact,self.Emax,self.del_max,self.R0, E0=self.E0, s=self.s,
            flag_costheta_delta_scale=self.flag_costheta_delta_scale, flag_costheta_Emax_shift=self.flag_costheta_Emax_shift)
        flag_elast=(rand(len(ref_frac))<ref_frac)
        flag_truesec=~(flag_elast)
        nel_emit=nel_impact*yiel

        return nel_emit, flag_elast, flag_truesec

    def impacts_on_surface(self, mass, nel_impact, x_impact, y_impact, z_impact, 
                                vx_impact, vy_impact, vz_impact, Norm_x, Norm_y, i_found,
                                v_impact_n, E_impact_eV, costheta_impact, nel_mp_th, flag_seg):

        
        nel_emit_tot_events, flag_elast, flag_truesec = self.SEY_process(nel_impact,E_impact_eV, costheta_impact, i_found)

        nel_replace = nel_emit_tot_events
        x_replace = x_impact.copy()
        y_replace = y_impact.copy()
        z_replace = z_impact.copy()
        vx_replace = vx_impact.copy()
        vy_replace = vy_impact.copy()
        vz_replace = vz_impact.copy()
        i_seg_replace = i_found.copy()

        # Handle elastics
        vx_replace[flag_elast], vy_replace[flag_elast] =  ee.specular_velocity(
                                                            vx_impact[flag_elast], vy_impact[flag_elast], 
                                                            Norm_x[flag_elast], Norm_y[flag_elast], v_impact_n[flag_elast]
                                                        )

        # true secondary
        N_true_sec = np.sum(flag_truesec)
        n_add_total = 0
        if N_true_sec > 0:

            n_add = np.zeros_like(flag_truesec, dtype=int)
            n_add[flag_truesec]=np.ceil(nel_replace[flag_truesec]/nel_mp_th)-1
            n_add[n_add<0]=0. #in case of underflow
            nel_replace[flag_truesec]=nel_replace[flag_truesec]/(n_add[flag_truesec]+1.)

            n_add_total = np.sum(n_add)

            # MPs to be replaced
            En_truesec_eV = ee.sec_energy_hilleret_model2(
                self.switch_no_increase_energy, N_true_sec, self.sigmafit, self.mufit, 
                self.E_th, E_impact_eV[flag_truesec], self.thresh_low_energy)

            vx_replace[flag_truesec], vy_replace[flag_truesec], vz_replace[flag_truesec] = self.angle_dist_func(
                N_true_sec, En_truesec_eV, Norm_x[flag_truesec], Norm_y[flag_truesec], mass)

            # Add new MPs
            if n_add_total != 0:
                # Clone MPs
                x_new_MPs = np.repeat(x_impact, n_add)
                y_new_MPs = np.repeat(y_impact, n_add)
                z_new_MPs = np.repeat(z_impact, n_add)
                norm_x_add = np.repeat(Norm_x, n_add)
                norm_y_add = np.repeat(Norm_y, n_add)
                nel_new_MPs = np.repeat(nel_replace, n_add)
                E_impact_eV_add = np.repeat(E_impact_eV, n_add)

                # Generate new MP properties, angles and energies
                En_truesec_eV_add = ee.sec_energy_hilleret_model2(
                    self.switch_no_increase_energy, n_add_total, self.sigmafit, self.mufit, 
                    self.E_th, E_impact_eV_add, self.thresh_low_energy)
                
                vx_new_MPs, vy_new_MPs, vz_new_MPs = self.angle_dist_func(
                    n_add_total, En_truesec_eV_add, norm_x_add, norm_y_add, mass)

                if flag_seg:
                    i_seg_new_MPs = np.repeat(i_found, n_add)
                else:
                    i_seg_new_MPs = None
        
        if n_add_total==0:
            nel_new_MPs = np.array([])
            x_new_MPs = np.array([])
            y_new_MPs = np.array([])
            z_new_MPs = np.array([])
            vx_new_MPs = np.array([])
            vy_new_MPs = np.array([])
            vz_new_MPs = np.array([])
            i_seg_new_MPs = np.array([])



        event_type = flag_truesec
        event_info = {}

        return nel_emit_tot_events, event_type, event_info,\
               nel_replace, x_replace, y_replace, z_replace, vx_replace, vy_replace, vz_replace, i_seg_replace,\
               nel_new_MPs, x_new_MPs, y_new_MPs, z_new_MPs, vx_new_MPs, vy_new_MPs, vz_new_MPs, i_seg_new_MPs



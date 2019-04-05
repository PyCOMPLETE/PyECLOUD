#-Begin-preamble-------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 7.7.0
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
#                           Eric Wulff
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

from numpy import sqrt, exp, take
from numpy.random import rand
import numpy as np
from sec_emission_model_ECLOUD import SEY_model_ECLOUD
from scipy.constants import e as qe

def yield_fun2(E, costheta, Emax, del_max, R0, E0):

    s = 1.35

    del_max_tilde = del_max * exp(0.5 * (1. - costheta))
    E_max_tilde = Emax * (1. + 0.7 * (1. - costheta))

    x = E / E_max_tilde

    true_sec = del_max_tilde * (s * x) / (s - 1. + x**s)
    reflected = R0 * ((sqrt(E) - sqrt(E + E0)) / (sqrt(E) + sqrt(E + E0)))**2.

    delta = true_sec + reflected
    ref_frac = 0. * delta
    mask_non_zero = (delta > 0)
    ref_frac[mask_non_zero] = reflected[mask_non_zero] / delta[mask_non_zero]

    return delta, ref_frac


class SEY_model_ECLOUD_non_unif(SEY_model_ECLOUD):
    
    def __init__(self, chamb, Emax, del_max, R0, E0=150.,
                 E_th=None, sigmafit=None, mufit=None,
                 switch_no_increase_energy=0, thresh_low_energy=None, secondary_angle_distribution=None,
                    ):

            if chamb.chamb_type != 'polyg':
                raise ValueError("""ECLOUD_nunif can be used only with chamb_type='polyg'!!!""")

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

            self.del_max_segments = np.float_(chamb.del_max_segments)
            self.R0_segments = np.float_(chamb.R0_segments)
            self.Emax_segments = np.float_(chamb.Emax_segments)

            self.del_max_segments[chamb.del_max_segments < 0.] = del_max
            self.R0_segments[chamb.R0_segments < 0.] = R0
            self.Emax_segments[chamb.Emax_segments < 0.] = Emax

            self.E0 = E0

            print 'Secondary emission model: ECLOUD non uniform E0=%f'%self.E0

    def SEY_process(self, nel_impact, E_impact_eV, costheta_impact, i_impact):

            Emax_mp = take(self.Emax_segments, i_impact)
            del_max_mp = take(self.del_max_segments, i_impact)
            R0_mp = take(self.R0_segments, i_impact)

            yiel, ref_frac = yield_fun2(E_impact_eV, costheta_impact, Emax_mp, del_max_mp, R0_mp, E0=self.E0)
            flag_elast = (rand(len(ref_frac)) < ref_frac)
            flag_truesec = ~(flag_elast)
            nel_emit = nel_impact * yiel

            return nel_emit, flag_elast, flag_truesec


class SEY_model_ECLOUD_non_unif_charging(SEY_model_ECLOUD_non_unif):
    
    def __init__(self, chamb, Emax, del_max, R0, E0=150.,
                 E_th=None, sigmafit=None, mufit=None,
                 switch_no_increase_energy=0, thresh_low_energy=None, secondary_angle_distribution=None,   
                 ):

        super(SEY_model_ECLOUD_non_unif_charging, self).__init__(chamb, Emax, del_max, R0, E0,
                    E_th, sigmafit, mufit,
                    switch_no_increase_energy, thresh_low_energy, secondary_angle_distribution,   
                    )
        print 'Secondary emission model: ECLOUD non uniform E0=%f, with charging'%self.E0
        
        self.chamb = chamb
        self.Q_segments = 0. * self.del_max_segments
        self.flag_charging =  np.int_(chamb.flag_charging)>0
        self.Q_max_segments = np.float_(chamb.Q_max_segments)
        self.EQ_segments = np.float_(chamb.EQ_segments)
        self.tau_segments = np.float_(chamb.tau_segments)


    def SEY_process(self, nel_impact, E_impact_eV, costheta_impact, i_impact):
        
        Emax_mp = take(self.Emax_segments, i_impact)
        del_max_mp = take(self.del_max_segments, i_impact)
        R0_mp = take(self.R0_segments, i_impact)

        yiel, ref_frac = yield_fun2(E_impact_eV, costheta_impact, Emax_mp, del_max_mp, R0_mp, E0=self.E0)
                 
        mask_charging = np.take(self.flag_charging, i_impact) 
        if np.any(mask_charging):
            Q_charging = np.take(self.Q_segments, i_impact[mask_charging])
            Q_max = np.take(self.Q_max_segments, i_impact[mask_charging])
            EQ = np.take(self.EQ_segments, i_impact[mask_charging])
            
            Q_charging[Q_charging<0.] = 0.
            Q_charging[Q_charging>Q_max] = Q_max[Q_charging>Q_max]
  
            yiel[mask_charging] = yiel[mask_charging] * (1. - Q_charging/Q_max) + (1. - np.exp(-E_impact_eV[mask_charging]/EQ))*(Q_charging/Q_max)
            
            # This would preserve also the ener
            # ref_frac[mask_charging] = ref_frac[mask_charging] * (1. - Q_charging/Q_max) + Q_charging/Q_max
        
        flag_elast = (rand(len(ref_frac)) < ref_frac)
        flag_truesec = ~(flag_elast)
   
        nel_emit = nel_impact * yiel

        for i_edg, flag_Q in enumerate(self.flag_charging):

            if flag_Q:
                mask_impact_here = (i_impact == i_edg)
                n_impact_here = np.sum(nel_impact[mask_impact_here])
                n_emit_here = np.sum(nel_emit[mask_impact_here])

                self.Q_segments[i_edg] += (n_impact_here - n_emit_here)*(-qe)/self.chamb.L_edg[i_edg]

        return nel_emit, flag_elast, flag_truesec

    def SEY_model_evol(self, Dt): 
        mask_evol = np.logical_and(self.flag_charging, self.tau_segments>0.)
        self.Q_segments[mask_evol]*=np.exp(-Dt/self.tau_segments[mask_evol])







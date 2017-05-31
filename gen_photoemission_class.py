#----------------------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#		           PyECLOUD Version 6.2.0
#
#
#     Author and contact:   Giovanni IADAROLA
#                           BE-ABP Group
#                           CERN
#                           CH-1211 GENEVA 23
#                           SWITZERLAND
#                           giovanni.iadarola@cern.ch
#
#                contact:   Giovanni RUMOLO
#                           BE-ABP Group
#                           CERN
#                           CH-1211 GENEVA 23
#                           SWITZERLAND
#                           giovanni.rumolo@cern.ch
#
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
#----------------------------------------------------------------------

from __future__ import division
import numpy.random as random
from numpy import floor, interp, pi, sin, cos, zeros, sqrt, sum

import scipy.io as sio
from scipy.constants import e as qe, m_e as me, c
import numpy as np

qm=qe/me

class photoemission:

    def __init__(self, inv_CDF_refl_photoem_file, k_pe_st, refl_frac, e_pe_sigma, e_pe_max,alimit, \
                x0_refl, y0_refl, out_radius, chamb, resc_fac, energy_distribution):

        print 'Start photoemission init.'

        if inv_CDF_refl_photoem_file == 'unif_no_file':
            self.flag_unif = True
        else:
            self.flag_unif = False
            dict_psi_inv_CDF = sio.loadmat(inv_CDF_refl_photoem_file)
            self.inv_CDF_refl = np.squeeze(dict_psi_inv_CDF['inv_CDF'].real)
            self.u_sam_CDF_refl = np.squeeze(dict_psi_inv_CDF['u_sam'].real)

        self.k_pe_st = k_pe_st
        self.refl_frac = refl_frac
        self.e_pe_sigma = e_pe_sigma
        self.e_pe_max = e_pe_max
        self.alimit = alimit
        self.x0_refl = x0_refl
        self.y0_refl = y0_refl
        self.out_radius = out_radius
        self.chamb = chamb
        self.resc_fac = resc_fac

        if y0_refl!=0.:
            raise ValueError('The case y0_refl!=0 is NOT IMPLEMETED yet!!!!')

        x0_refl_np_arr = np.array([x0_refl])
        y0_refl_np_arr = np.array([y0_refl])
        if np.any(self.chamb.is_outside(x0_refl_np_arr, y0_refl_np_arr)):
            raise ValueError('x0_refl, y0_refl is outside of the chamber!')

        # Convert normal distribution parameters to lognormal distribution parameters
        if energy_distribution == 'lognormal':
            e_pe_sigma_logn = np.sqrt(np.log(e_pe_sigma**2/e_pe_max**2 + 1.))
            e_pe_max_logn  = np.log(e_pe_max) - e_pe_sigma_logn**2/2.
            def get_energy(N_int_new_MP):
                return random.lognormal(e_pe_max_logn, e_pe_sigma_logn, N_int_new_MP)

        elif energy_distribution == 'gaussian':
            def get_energy(N_int_new_MP):
                En_gen = random.randn(N_int_new_MP)*e_pe_sigma + e_pe_max
                flag_negat=(En_gen<0.)
                N_neg=sum(flag_negat);
                while(N_neg>0):
                    En_gen[flag_negat]=random.randn(N_neg)*e_pe_sigma + e_pe_max   #in eV
                    flag_negat=(En_gen<0.)
                    N_neg=sum(flag_negat);
                return En_gen
        elif energy_distribution == 'rect':
            def get_energy(N_int_new_MP):
                return e_pe_max + (random.rand(N_int_new_MP)-0.5)*e_pe_sigma
        elif energy_distribution == 'mono':
            def get_energy(N_int_new_MP):
                return np.ones(N_int_new_MP)*e_pe_max
        else:
            raise ValueError('Energy distribution not specified')
        self.get_energy = get_energy
        print 'Done photoemission init. Energy distribution: %s' % energy_distribution

    def generate(self, MP_e, lambda_t, Dt):

        k_pe=self.k_pe_st*c
        #determine the number of MPs to be generated
        DNel=k_pe*lambda_t*Dt

        N_new_MP=DNel/MP_e.nel_mp_ref
        Nint_new_MP=floor(N_new_MP)
        rest=N_new_MP-Nint_new_MP
        Nint_new_MP=Nint_new_MP+int(random.rand()<rest)
        Nint_new_MP=int(Nint_new_MP)


        if Nint_new_MP>0:
            #generate appo x_in and x_out
            x_in = zeros(Nint_new_MP)
            y_in = zeros(Nint_new_MP)
            x_out = zeros(Nint_new_MP)
            y_out = zeros(Nint_new_MP)

            #for each one generate flag refl
            refl_flag=(random.rand(Nint_new_MP)<self.refl_frac)
            gauss_flag=~refl_flag

            #generate psi for refl. photons generation
            N_refl=sum(refl_flag)
            if N_refl>0:
                u_gen=random.rand(N_refl,1);
                if self.flag_unif:
                    psi_gen = 2.*pi*u_gen
                    x_out[refl_flag]=self.out_radius*cos(psi_gen)
                    y_out[refl_flag]=self.out_radius*sin(psi_gen)
                else:
                    psi_gen=interp(u_gen,self.u_sam_CDF_refl, self.inv_CDF_refl);
                    x_in[refl_flag]=self.x0_refl
                    x_out[refl_flag]=-2.*self.out_radius*cos(psi_gen)+self.x0_refl
                    y_out[refl_flag]=2.*self.out_radius*sin(psi_gen)


            #generate theta for nonreflected photon generation
            N_gauss=sum(gauss_flag)
            if N_gauss>0:
                #theta_gen=self.alimit*randn(N_gauss)
                theta_gen = random.normal(0, self.alimit, N_gauss)
                x_out[gauss_flag]=self.out_radius*cos(theta_gen)
                y_out[gauss_flag]=self.out_radius*sin(theta_gen)

            #generate points and normals
            x_int, y_int, _, Norm_x, Norm_y, i_found= self.chamb.impact_point_and_normal(x_in, y_in, 0*x_in,
                                                                                          x_out, y_out, 0*x_out, resc_fac=self.resc_fac)


            #generate energies (the same distr. for all photoelectr.)
            En_gen = self.get_energy(N_new_MP) #in eV


            # generate velocities like in impact managment
            v_gen_mod=sqrt(2.*qm*En_gen);

            sin_theta_p=random.rand(Nint_new_MP);
            cos_theta_p=sqrt(1.-sin_theta_p*sin_theta_p);
            phi_p=random.rand(Nint_new_MP)*2.*pi;
            sin_phi_p=sin(phi_p);
            cos_phi_p=cos(phi_p);

            vx_gen=v_gen_mod*\
                (cos_theta_p*Norm_x+sin_theta_p*sin_phi_p*Norm_y);
            vy_gen=v_gen_mod*\
                (cos_theta_p*Norm_y-sin_theta_p*sin_phi_p*Norm_x);
            vz_gen=v_gen_mod*(sin_theta_p*cos_phi_p);


            MP_e.x_mp[MP_e.N_mp:MP_e.N_mp+Nint_new_MP]=x_int;#Be careful to the indexing when translating to python
            MP_e.y_mp[MP_e.N_mp:MP_e.N_mp+Nint_new_MP]=y_int;
            MP_e.z_mp[MP_e.N_mp:MP_e.N_mp+Nint_new_MP]=0.;#randn(Nint_new_MP,1);
            MP_e.vx_mp[MP_e.N_mp:MP_e.N_mp+Nint_new_MP]=vx_gen
            MP_e.vy_mp[MP_e.N_mp:MP_e.N_mp+Nint_new_MP]=vy_gen
            MP_e.vz_mp[MP_e.N_mp:MP_e.N_mp+Nint_new_MP]=vz_gen
            MP_e.nel_mp[MP_e.N_mp:MP_e.N_mp+Nint_new_MP]=MP_e.nel_mp_ref;

            MP_e.N_mp=int(MP_e.N_mp+Nint_new_MP);


        return MP_e


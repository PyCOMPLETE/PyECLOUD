#----------------------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 6.4.0
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

from __future__ import division, print_function
import numpy as np
import numpy.random as random

import scipy.io as sio
import scipy.stats as stats
from scipy.constants import e as qe, m_e as me, c

import sec_emission

qm=qe/me

class photoemission(object):

    def __init__(self, inv_CDF_refl_photoem_file, k_pe_st, refl_frac, e_pe_sigma, e_pe_max,alimit, x0_refl,
                 y0_refl, out_radius, chamb, resc_fac, energy_distribution, photoelectron_angle_distribution):

        print('Start photoemission init.')

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
        self.angle_dist_func = sec_emission.get_angle_dist_func(photoelectron_angle_distribution)

        if y0_refl!=0.:
            raise ValueError('The case y0_refl!=0 is NOT IMPLEMETED yet!!!!')

        x0_refl_np_arr = np.array([x0_refl])
        y0_refl_np_arr = np.array([y0_refl])
        if np.any(self.chamb.is_outside(x0_refl_np_arr, y0_refl_np_arr)):
            raise ValueError('x0_refl, y0_refl is outside of the chamber!')

        self.get_energy = self.get_energy_distribution_func(energy_distribution, e_pe_sigma, e_pe_max)
        print('Done photoemission init. Energy distribution: %s' % energy_distribution)

    def get_energy_distribution_func(self, energy_distribution, e_pe_sigma, e_pe_max):
        if energy_distribution == 'lognormal':

            def get_energy(N_int_new_MP):
                return random.lognormal(e_pe_max, e_pe_sigma, N_int_new_MP)

        elif energy_distribution == 'gaussian':

            def get_energy(N_int_new_MP):
                En_gen = random.randn(N_int_new_MP)*e_pe_sigma + e_pe_max
                flag_negat = (En_gen<0.)
                N_neg = sum(flag_negat)
                while(N_neg>0):
                    En_gen[flag_negat] = random.randn(N_neg)*e_pe_sigma + e_pe_max   #in eV
                    flag_negat = (En_gen<0.)
                    N_neg = sum(flag_negat)
                return En_gen

        elif energy_distribution == 'rect':

            def get_energy(N_int_new_MP):
                return e_pe_max + (random.rand(N_int_new_MP)-0.5)*e_pe_sigma

        elif energy_distribution == 'mono':

            def get_energy(N_int_new_MP):
                return np.ones(N_int_new_MP)*e_pe_max

        elif energy_distribution == 'lorentz':

            xx_min = stats.cauchy.cdf(0, e_pe_max, e_pe_sigma)
            xx_max = 1

            def get_energy(N_int_new_MP):
                xx_rand = random.rand(N_int_new_MP)*(xx_max-xx_min) + xx_min
                return stats.cauchy.ppf(xx_rand, e_pe_max, e_pe_sigma)
        else:
            raise ValueError('Energy distribution not specified')

        return get_energy

    def get_number_new_mps(self, lambda_t, Dt, nel_mp_ref):
        k_pe = self.k_pe_st*c
        DNel = k_pe*lambda_t*Dt

        N_new_MP=DNel/nel_mp_ref
        rest, Nint_new_MP = np.modf(N_new_MP)
        return int(Nint_new_MP+int(random.rand()<rest))

    def gen_energy_and_set_MPs(self, Nint_new_MP, x_in, y_in, x_out, y_out, MP_e):
        #generate points and normals

        z_in = z_out = np.zeros_like(x_out)

        x_int, y_int, _, Norm_x, Norm_y, i_found= self.chamb.impact_point_and_normal(
            x_in, y_in, z_in,x_out, y_out, z_out, resc_fac=self.resc_fac)

        #generate energies (the same distr. for all photoelectr.)
        En_gen = self.get_energy(Nint_new_MP) #in eV
        # generate velocities like in impact managment
        vx_gen, vy_gen, vz_gen = self.angle_dist_func(Nint_new_MP, En_gen, Norm_x, Norm_y)

        MP_e.add_new_MPs(Nint_new_MP, MP_e.nel_mp_ref, x_int, y_int, 0., vx_gen, vy_gen, vz_gen)

    def generate(self, MP_e, lambda_t, Dt):

        Nint_new_MP = self.get_number_new_mps(lambda_t, Dt, MP_e.nel_mp_ref)

        if Nint_new_MP>0:
            #generate appo x_in and x_out
            x_in = np.zeros(Nint_new_MP)
            y_in = np.zeros(Nint_new_MP)
            x_out = np.zeros(Nint_new_MP)
            y_out = np.zeros(Nint_new_MP)

            #for each one generate flag refl
            refl_flag=(random.rand(Nint_new_MP)<self.refl_frac)
            gauss_flag=~refl_flag

            #generate psi for refl. photons generation
            N_refl=np.sum(refl_flag)
            if N_refl>0:
                u_gen=random.rand(N_refl)
                if self.flag_unif:
                    psi_gen = 2.*np.pi*u_gen
                    x_out[refl_flag]=self.out_radius*np.cos(psi_gen)
                    y_out[refl_flag]=self.out_radius*np.sin(psi_gen)
                else:
                    psi_gen = np.interp(u_gen, self.u_sam_CDF_refl, self.inv_CDF_refl)
                    x_in[refl_flag]=self.x0_refl
                    x_out[refl_flag]=-2.*self.out_radius*np.cos(psi_gen)+self.x0_refl
                    y_out[refl_flag]=2.*self.out_radius*np.sin(psi_gen)

            #generate theta for nonreflected photon generation
            N_gauss=np.sum(gauss_flag)
            if N_gauss>0:
                #theta_gen=self.alimit*randn(N_gauss)
                theta_gen = random.normal(0, self.alimit, N_gauss)
                x_out[gauss_flag]=self.out_radius*np.cos(theta_gen)
                y_out[gauss_flag]=self.out_radius*np.sin(theta_gen)

            self.gen_energy_and_set_MPs(Nint_new_MP, x_in, y_in, x_out, y_out, MP_e)

        return MP_e


class photoemission_from_file(photoemission):
    def __init__(self, photoemission_file, chamb, resc_frac, energy_distribution, e_pe_sigma, e_pe_max, k_pe_st, out_radius, photoelectron_angle_distribution):
        print('Start photoemission init from file %s.' % photoemission_file)

        mat = sio.loadmat(photoemission_file)
        self.inv_CDF = mat['inv_CDF']
        self.angles = mat['angles']
        self.k_pe_st = k_pe_st
        self.out_radius = out_radius
        self.chamb = chamb
        self,resc_frac = resc_frac

        self.get_energy = self.get_energy_distribution_func(energy_distribution, e_pe_sigma, e_pe_max)
        self.angle_dist_func = sec_emission.get_angle_dist_func(photoelectron_angle_distribution)
        print('Done photoemission init')

    def generate(self, MP_e, lambda_t, Dt):

        Nint_new_MP = self.get_number_new_mps(lambda_t, Dt, MP_e.nel_mp_ref)
        if Nint_new_MP>0:

            cdf_gen = random.rand(Nint_new_MP)
            theta_gen = np.interp(cdf_gen, self.inv_CDF, self.angles)

            x_out = self.out_radius*np.cos(theta_gen)
            y_out = self.out_radius*np.sin(theta_gen)

            x_in = y_in = np.zeros(Nint_new_MP)
            self.gen_energy_and_set_MPs(Nint_new_MP, x_in, y_in, x_out, y_out, MP_e)

        return MP_e


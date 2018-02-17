#-Begin-preamble-------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 6.7.1
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

from __future__ import division
import time
import numpy as np
import numpy.random as random
from scipy.constants import e as qe, m_e
import scipy.stats as stats

qm = qe/m_e

# Secondary electron secondaries
def sec_energy_hilleret_model2(switch_no_increase_energy, Ngen, sigmafit, mufit, E_th, En_impact_eV, thresh_low_energy):

    if switch_no_increase_energy==0:
        en_eV=random.lognormal(mufit,sigmafit,Ngen)
        flag_above_th=(en_eV>E_th)
        Nabove_th = np.sum(flag_above_th)

        while Nabove_th>0:
            en_eV[flag_above_th]=random.lognormal(mufit,sigmafit,Nabove_th)

            flag_above_th=(en_eV>E_th)
            Nabove_th = np.sum(flag_above_th)

    elif switch_no_increase_energy==1:

        en_eV = np.zeros_like(En_impact_eV, dtype=float)

        flag_low_energy = En_impact_eV<thresh_low_energy
        flag_high_energy = ~(flag_low_energy)
        N_low_ene = np.sum(flag_low_energy)
        N_high_ene = np.sum(flag_high_energy)

        #generate low energy
        en_eV_le=random.randn(N_low_ene)   #in eV
        flag_negat=np.logical_or(en_eV_le<0., en_eV_le>4.)
        N_neg=np.sum(flag_negat)
        while(N_neg>0):
            en_eV_le[flag_negat]=random.randn(N_neg)  #in eV
            flag_negat=np.logical_or(en_eV_le<0., en_eV_le>4.)
            N_neg=np.sum(flag_negat)
        sigma_le=En_impact_eV[flag_low_energy]/4.
        en_eV_le=(en_eV_le+2.)*sigma_le

        #generate high energy
        en_eV_he=random.lognormal(mufit,sigmafit,N_high_ene)

        flag_above_th = np.logical_or(en_eV_he>E_th,(en_eV_he-En_impact_eV[flag_high_energy])>0)
        Nabove_th = np.sum(flag_above_th)

        while Nabove_th>0:
            en_eV_he[flag_above_th]=random.lognormal(mufit,sigmafit,Nabove_th)
            flag_above_th = np.logical_or(en_eV_he>E_th,(en_eV_he-En_impact_eV[flag_high_energy])>0)
            Nabove_th = np.sum(flag_above_th)

        en_eV[flag_high_energy]=en_eV_he
        en_eV[flag_low_energy]=en_eV_le

    return en_eV


# Angle of generated electrons

# This function correctly implements a 3D cosine distribution.
# The sin(theta) factor that occurs in 3D spherical integrals is taken into account.
# See also the paper from J.Greenwood,
# "The correct and incorrect generation of a cosine distribution of scattered
# particles for Monte-Carlo modelling of vacuum systems"
# http://www.sciencedirect.com/science/article/pii/S0042207X02001732
#
# And further in the proceedings of ECLOUD '02 workshop p.105:
# Rumolo and Zimmermann: "Electron-Cloud Simulations: Build Up and related effects.
# https://cds.cern.ch/record/537336?ln=en

# Fixed behavior
def velocities_angle_cosine_3D(N_new_MP, En_gen, Norm_x, Norm_y):
    sin_theta_p = np.sqrt(random.rand(N_new_MP))
    return _velocities_angle(N_new_MP, En_gen, Norm_x, Norm_y, sin_theta_p)

# This has been the behavior of the code until the error was spotted.
def velocities_angle_cosine_2D(N_new_MP, En_gen, Norm_x, Norm_y):
    sin_theta_p = random.rand(N_new_MP)
    return _velocities_angle(N_new_MP, En_gen, Norm_x, Norm_y, sin_theta_p)

# Avoid code duplication
def _velocities_angle(N_new_MP, En_gen, Norm_x, Norm_y, sin_theta_p):
    v_gen_mod=np.sqrt(2.*qm*En_gen)

    phi_p = random.rand(N_new_MP)*2*np.pi
    sin_phi_p = np.sin(phi_p)
    cos_phi_p = np.cos(phi_p)

    cos_theta_p = np.sqrt(1-sin_theta_p**2)

    vx_gen = v_gen_mod*(cos_theta_p*Norm_x+sin_theta_p*sin_phi_p*Norm_y)
    vy_gen = v_gen_mod*(cos_theta_p*Norm_y-sin_theta_p*sin_phi_p*Norm_x)
    vz_gen = v_gen_mod*(sin_theta_p*cos_phi_p)

    return vx_gen, vy_gen, vz_gen

# Interface
def get_angle_dist_func(string):
    if string == 'cosine_3D':
        print('Using cosine_3D emission angle distribution.')
        return velocities_angle_cosine_3D
    elif string == 'cosine_2D':
        print("""
Warning! The 2D emission angle distribution is used!
The 'cosine_3D' distribution is more appropriate. This can be enabled by
setting in the machine parameter input file:
"photoelectron_angle_distribution = 'cosine_3D'"
and in the secondary emission input file:
"secondary_angle_distribution = 'cosine_3D'"
For more info, see presentation by P. Dijkstal on the angle of emission 
of generated electrons (https://indico.cern.ch/event/673160/).
""")
        time.sleep(3)
        return velocities_angle_cosine_2D
    else:
        raise ValueError("""
The emission angle distribution must be specified!
To use the cosine_3D distribution (most appropriate) set in the
machine parameter input file:
"photoelectron_angle_distribution = 'cosine_3D'"
and in the secondary emission input file:
"secondary_angle_distribution = 'cosine_3D'"
In case you want to use the cosine_2D distribution, to have the same
behavior as in PyECLOUD 6.6.0 or earlier replace 'cosine_3D' with
'cosine_2D'.
For more info, see presentation by P. Dijkstal on the angle of emission
of generated electrons (https://indico.cern.ch/event/673160/).
""")


# Photoelectron energy classes
class _gen_energy_base(object):
    def __init__(self, e_pe_sigma, e_pe_max):
        self.e_pe_sigma = e_pe_sigma
        self.e_pe_max = e_pe_max

class _lognormal(_gen_energy_base):
    def __call__(self, N_int_new_MP):
        return random.lognormal(self.e_pe_max, self.e_pe_sigma, N_int_new_MP)

class _gaussian(_gen_energy_base):
    def __call__(self, N_int_new_MP):
        En_gen = random.randn(N_int_new_MP)*self.e_pe_sigma + self.e_pe_max
        flag_negat = (En_gen<0.)
        N_neg = np.sum(flag_negat)
        while(N_neg>0):
            En_gen[flag_negat] = random.randn(N_neg)*self.e_pe_sigma +self.e_pe_max   #in eV
            flag_negat = (En_gen<0.)
            N_neg = np.sum(flag_negat)
        return En_gen

class _rect(_gen_energy_base):
    def __call__(self, N_int_new_MP):
        return self.e_pe_max + (random.rand(N_int_new_MP)-0.5)*self.e_pe_sigma

class _mono(_gen_energy_base):
    def __call__(self, N_int_new_MP):
        return np.ones(N_int_new_MP)*self.e_pe_max

class _lorentz(_gen_energy_base):
    def __init__(self, e_pe_sigma, e_pe_max):
        self.e_pe_sigma = e_pe_sigma
        self.e_pe_max = e_pe_max
        self.xx_min = stats.cauchy.cdf(0, e_pe_max, e_pe_sigma)
        self.xx_max = 1 # set this to something else if you want to cut

    def __call__(self, N_int_new_MP):
        xx_rand = random.rand(N_int_new_MP)*(self.xx_max-self.xx_min) + self.xx_min
        return stats.cauchy.ppf(xx_rand, self.e_pe_max, self.e_pe_sigma)

# Interface
def get_energy_distribution_func(energy_distribution, e_pe_sigma, e_pe_max):

    if energy_distribution == 'lognormal':
        get_energy = _lognormal
    elif energy_distribution == 'gaussian':
        get_energy = _gaussian
    elif energy_distribution == 'rect':
        get_energy = _rect
    elif energy_distribution == 'mono':
        get_energy = _mono
    elif energy_distribution == 'lorentz':
        get_energy = _lorentz
    else:
        raise ValueError('Energy distribution %s is invalid!' % energy_distribution)

    return get_energy(e_pe_sigma, e_pe_max)



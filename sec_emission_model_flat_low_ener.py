#-Begin-preamble-------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 8.5.0
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
#                           Lorenzo Giacomel
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

from numpy import sqrt, exp, cos, pi, logical_and
from numpy.random import rand
from .sec_emission_model_ECLOUD import SEY_model_ECLOUD


def yield_fun2(E, costheta, Emax, del_max, R0):

    s = 1.35

    del_max_tilde = del_max * exp(0.5 * (1. - costheta))
    E_max_tilde = Emax * (1. + 0.7 * (1. - costheta))

    x = E / E_max_tilde

    true_sec = del_max_tilde * (s * x) / (s - 1. + x**s)

    delta = true_sec

    mask_change = logical_and(delta < R0, E < Emax)
    delta[mask_change] = R0
    reflected = delta - true_sec

    ref_frac = reflected / delta

    return delta, ref_frac


class SEY_model_flat_le(SEY_model_ECLOUD):
    def __init__(self, Emax, del_max, R0,
                 E_th=None, sigmafit=None, mufit=None,
                 switch_no_increase_energy=0, thresh_low_energy=None, secondary_angle_distribution=None,
                    ):

            self.E_th = E_th
            self.sigmafit = sigmafit
            self.mufit = mufit
            self.switch_no_increase_energy = switch_no_increase_energy
            self.thresh_low_energy = thresh_low_energy
            self.secondary_angle_distribution = secondary_angle_distribution

            if secondary_angle_distribution is not None:
                from . import electron_emission
                self.angle_dist_func = electron_emission.get_angle_dist_func(secondary_angle_distribution)
            else:
                self.angle_dist_func = None

            self.Emax = Emax
            self.del_max = del_max
            self.R0 = R0
            print('Secondary emission model: Flat Low Energy ')

    def SEY_process(self, nel_impact, E_impact_eV, costheta_impact, i_impact):
            yiel, ref_frac = yield_fun2(E_impact_eV, costheta_impact, self.Emax, self.del_max, self.R0)
            flag_elast = (rand(len(ref_frac)) < ref_frac)
            flag_truesec = ~(flag_elast)
            nel_emit = nel_impact * yiel

            return nel_emit, flag_elast, flag_truesec


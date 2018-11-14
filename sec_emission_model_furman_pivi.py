#-Begin-preamble-------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 7.6.0
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
import numpy.random as random
from sec_emission_model_ECLOUD import SEY_model_ECLOUD


def energy_rediffused(self, E0):
    randn = random.rand(len(E0))
    return randn**(1 / (self.q + 1)) * E0


def energy_trueSecondary(self, E0):
    u = random.rand(len(E0))
    E_out = self.eHat0
    return E_out


class SEY_model_furman_pivi(SEY_model_ECLOUD):
    def __init__(self, Emax, del_max, R0,
                 E_th=None, sigmafit=None, mufit=None,
                 switch_no_increase_energy=0, thresh_low_energy=None, secondary_angle_distribution=None,
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

        print 'Secondary emission model: Furman-Pivi E0=%.4f s=%.4f' % (self.E0, self.s)

    def SEY_process(self, nel_impact, E_impact_eV, costheta_impact, i_impact):
        # Furman-Pivi algorithm
        # (1): Compute emission angles and energy

        # Already implemented in the impact_man class.

        # (2): Compute delta_e, delta_r, delta_ts
        delta_e, delta_r, delta_ts = self._yield_fun_furman_pivi(self, E_impact_eV, costheta_impact)

        # (3): Generate probability of number electrons created

        # Emission probability per penetrated electron

        # Decide on type
        rand = random.rand(E_impact_eV.size)
        flag_rediffused = rand < delta_r
        flag_backscattered = np.logical_and(~flag_rediffused, rand < delta_r + delta_e)
        flag_truesec = np.logical_and(~flag_rediffused, ~flag_backscattered)

        # Reflected or backscattered electrons have yield 1 by definition.
        delta = np.ones_like(E_impact_eV, dtype=float)
        # True secondary part has to be adjusted accordingly.
        delta[flag_truesec] = delta_ts[flag_truesec] / (1. - delta_r[flag_truesec] - delta_e[flag_truesec])  # Eq. (39) in FP paper

        # (4): Generate number of secondaries for every impact
        # In impacts_on_surface

        # (5): Delete if n = 0
        # Done automatically by the MP system.

        # (6): Generate energy:
        # In impacts_on_surface

        nel_emit = delta * nel_impact

        return nel_emit, flag_backscattered, flag_rediffused, flag_truesec

    # def _yield_fun_furman_pivi(self, E_impact_eV, costheta_impact):
    #
    #     # elastic
    #     delta_e = self.p1EInf + (self.p1Ehat - self.p1EInf) * np.exp(-(np.abs(E_impact_eV - self.eEHat) / self.w)**self.p / self.p) * (1. + self.e1 * (1. - costheta_impact**self.e2))
    #     delta_r = self.p1RInf * (1. - np.exp(-(E_impact_eV / self.eR)**self.r)) * (1. + self.r1 * (1. - costheta_impact**self.r2))
    #     delta_ts = self.deltaTSHat * self._D(E_impact_eV / (self.eHat0 * (1. + self.t3 * (1. - np.cos(costheta_impact)**self.t4)))) * (1. + self.t1 * (1. - costheta_impact**self.t2))
    #
    #     return delta_e, delta_r, delta_ts

    def _D(self, x):
        s = self.s
        return s * x / (s - 1 + x**s)

    def _yield_fun_furman_pivi(self, E, costheta):
        delta_e = self._delta_e
        delta_r = self._delta_r
        delta_ts = self._delta_ts
        return delta_e, delta_r, delta_ts

    def _delta_e(self, E_impact_eV, costheta_impact):
        """
        Backscattered electrons (elastically scattered).
        (25) in FP paper.
        """

        exp_factor = -(np.abs(E_impact_eV - self.eEHat) / self.w)**self.p / self.p
        delta_e0 = self.p1EInf + (self.p1Ehat - self.p1EInf) * np.exp(exp_factor)
        angular_factor = 1. + self.e1 * (1. - costheta_impact**self.e2)

        return delta_e0 * angular_factor

    def _delta_r(self, E_impact_eV, costheta_impact):
        """
        Rediffused electrons (not in ECLOUD model).
        (28) in FP paper.
        """

        exp_factor = -(E_impact_eV / self.eR)**self.r
        delta_r0 = self.p1RInf * (1. - np.exp(exp_factor))
        angular_factor = 1. + self.r1 * (1. - costheta_impact**self.r2)

        return delta_r0 * angular_factor

    def _delta_ts(self, E_impact_eV, costheta_impact):
        """
        True secondaries.
        (31) in FP paper.
        """

        eHat = self.eHat0 * (1. + self.t3 * (1. - np.cos(costheta_impact)**self.t4))
        delta_ts0 = self.deltaTSHat * self._D(E_impact_eV / eHat)
        angular_factor = 1. + self.t1 * (1. - costheta_impact**self.t2)

        return delta_ts0 * angular_factor

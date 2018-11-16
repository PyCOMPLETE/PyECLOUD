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

import numpy as np
import numpy.random as random
import electron_emission as ee
# from scipy.special import gamma

# def inverse_CDF_ts_energy(delta_e, delta_r, delta_ts, n, pn, epsn):
#     delta_prime_ts = delta_ts / (1. - delta_e - delta_r)
#     PnTS_prime = delta_prime_ts**n / np.math.factorial(n) * np.exp(-delta_prime_ts)
#     Fn = PnTS_prime / ()(eps[n-1]**pn[n-1] * gamma(pn[n-1]))**n * )


class SEY_model_furman_pivi():
    def __init__(self,  # eHat0, deltaTSHat,
                 E_th=None, sigmafit=None, mufit=None,
                 switch_no_increase_energy=0, thresh_low_energy=None, secondary_angle_distribution=None,
                 s=1.54
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

        # self.eHat0 = eHat0
        # self.deltaTSHat = deltaTSHat
        # self.s = s

        print 'Secondary emission model: Furman-Pivi s=%.4f' % (self.s)

    def SEY_process(self, nel_impact, E_impact_eV, costheta_impact, i_impact):
        # Furman-Pivi algorithm
        # (1): Compute emission angles and energy

        # Already implemented in the impact_man class.

        # (2): Compute delta_e, delta_r, delta_ts
        delta_e, delta_r, delta_ts = self._yield_fun_furman_pivi(E_impact_eV, costheta_impact)

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

    def _yield_fun_furman_pivi(self, E, costheta):
        delta_e = self._delta_e(E, costheta)
        delta_r = self._delta_r(E, costheta)
        delta_ts = self._delta_ts(E, costheta)
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

        eHat = self.eHat0 * (1. + self.t3 * (1. - costheta_impact**self.t4))
        delta_ts0 = self.deltaTSHat * self._D(E_impact_eV / eHat)
        angular_factor = 1. + self.t1 * (1. - costheta_impact**self.t2)

        return delta_ts0 * angular_factor

    def _D(self, x):
        s = self.s
        return s * x / (s - 1 + x**s)

    def energy_rediffused(self, E0):
        uu = random.rand(len(E0))
        return uu**(1 / (self.q + 1)) * E0  # Inverse transform sampling of (29) in FP paper

    def energy_trueSecondary(self, E0):
        uu = random.rand(len(E0))

        E_out = self.eHat0
        return E_out

    def impacts_on_surface(self, mass, nel_impact, x_impact, y_impact, z_impact,
                           vx_impact, vy_impact, vz_impact, Norm_x, Norm_y, i_found,
                           v_impact_n, E_impact_eV, costheta_impact, nel_mp_th, flag_seg):

        nel_emit_tot_events, flag_backscattered, flag_rediffused, flag_truesec = self.SEY_process(nel_impact, E_impact_eV, costheta_impact, i_found)

        nel_replace = nel_emit_tot_events.copy()
        x_replace = x_impact.copy()
        y_replace = y_impact.copy()
        z_replace = z_impact.copy()
        vx_replace = vx_impact.copy()
        vy_replace = vy_impact.copy()
        vz_replace = vz_impact.copy()
        if i_found is not None:
            i_seg_replace = i_found.copy()
        else:
            i_seg_replace = i_found

        # Handle elastics
        vx_replace[flag_backscattered], vy_replace[flag_backscattered] = ee.specular_velocity(
            vx_impact[flag_backscattered], vy_impact[flag_backscattered],
            Norm_x[flag_backscattered], Norm_y[flag_backscattered], v_impact_n[flag_backscattered]
        )

        # Rediffused
        En_rediffused_eV = self.energy_rediffused(E_impact_eV[flag_rediffused])
        N_rediffused = np.sum(flag_rediffused)
        vx_replace[flag_rediffused], vy_replace[flag_rediffused], vz_replace[flag_rediffused] = self.angle_dist_func(
            N_rediffused, En_rediffused_eV, Norm_x[flag_rediffused], Norm_y[flag_rediffused], mass)

        # True secondary
        N_true_sec = np.sum(flag_truesec)
        n_add_total = 0
        if N_true_sec > 0:

            n_add = np.zeros_like(flag_truesec, dtype=int)
            n_add[flag_truesec] = np.ceil(nel_replace[flag_truesec] / nel_mp_th) - 1
            n_add[n_add < 0] = 0.  # in case of underflow
            nel_replace[flag_truesec] = nel_replace[flag_truesec] / (n_add[flag_truesec] + 1.)

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

        if n_add_total == 0:
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


class SEY_model_FP_Cu(SEY_model_furman_pivi):

    # Parameters for backscattered (elastically scattered) electrons
    # (25) in FP paper
    p1EInf = 0.02      # Minimum probability of elastic scattering (at infinite energy)
    p1Ehat = 0.496     # Peak probability
    eEHat = 0.        # Peak energy
    w = 60.86     # Exponential factor 1
    p = 1.        # Exponential factor 2
    # (47a)                 # Angular factors
    e1 = 0.26
    e2 = 2.
    # (26)
    sigmaE = 2.

    # Parameters for rediffused electrons
    # (28)
    p1RInf = 0.2       # Minimum probability of rediffused scattering (at infinite energy)
    eR = 0.041     # Peak energy
    r = 0.104     # Exponential factor
    # (29)
    q = 0.5
    # (47b)                 # Angular factors
    r1 = 0.26
    r2 = 2.

    # Parameters for true secondaries
    # (31)
    deltaTSHat = 1.8848    # Maximum probability of secondaries
    eHat0 = 276.8     # Peak enery
    # (32)
    s = 1.54      # Form factor of fitting curve
    # (48a)                 # Angular factors
    t1 = 0.66
    t2 = 0.8
    # (48b)
    t3 = 0.7
    t4 = 1.

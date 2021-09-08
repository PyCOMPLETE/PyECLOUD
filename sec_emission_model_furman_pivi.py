#-Begin-preamble-------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 8.5.1
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

import numpy as np
import numpy.random as random
import scipy
from scipy.special import gamma
from scipy.special import gammainc
from scipy.special import gammaincinv
from scipy.special import binom
from scipy.special import erf
from scipy.special import erfinv
from . import electron_emission as ee

_factorial = np.array([1,
                       1,
                       2,
                       6,
                       24,
                       120,
                       720,
                       5040,
                       40320,
                       362880,
                       3628800,
                       39916800,
                       479001600,
                       6227020800,
                       87178291200,
                       1307674368000,
                       20922789888000,
                       355687428096000,
                       6402373705728000,
                       121645100408832000,
                       2432902008176640000])


def factorial(n):
    return _factorial[np.array(n)]


class SEY_model_furman_pivi():

    event_types = {
        0: 'elast',
        1: 'true',
        2: 'rediff',
        3: 'absorb',
    }

    def __init__(self, furman_pivi_surface,
                 E_th=None, sigmafit=None, mufit=None,
                 switch_no_increase_energy=0, thresh_low_energy=None, secondary_angle_distribution=None,
                 flag_costheta_delta_scale=True, flag_costheta_Emax_shift=True):

        self.E_th = E_th
        self.sigmafit = sigmafit
        self.mufit = mufit
        self.switch_no_increase_energy = switch_no_increase_energy
        self.thresh_low_energy = thresh_low_energy
        self.secondary_angle_distribution = secondary_angle_distribution

        if secondary_angle_distribution is not None:
            self.angle_dist_func = ee.get_angle_dist_func(secondary_angle_distribution)
        else:
            self.angle_dist_func = None

        self.flag_costheta_delta_scale = flag_costheta_delta_scale
        self.flag_costheta_Emax_shift = flag_costheta_Emax_shift

        # General FP model parameters
        self.use_modified_sigmaE = furman_pivi_surface['use_modified_sigmaE']
        self.use_ECLOUD_theta0_dependence = furman_pivi_surface['use_ECLOUD_theta0_dependence']
        self.use_ECLOUD_energy = furman_pivi_surface['use_ECLOUD_energy']
        self.conserve_energy = furman_pivi_surface['conserve_energy']
        self.choice = furman_pivi_surface['choice']
        self.M_cut = furman_pivi_surface['M_cut']
        self.p_n = furman_pivi_surface['p_n']
        self.eps_n = furman_pivi_surface['eps_n']
        self.exclude_rediffused = furman_pivi_surface['exclude_rediffused']

        # Parameters for backscattered (elastically scattered) electrons
        self.p1EInf = furman_pivi_surface['p1EInf']
        self.p1Ehat = furman_pivi_surface['p1Ehat']
        self.eEHat = furman_pivi_surface['eEHat']
        self.w = furman_pivi_surface['w']
        self.p = furman_pivi_surface['p']
        self.e1 = furman_pivi_surface['e1']
        self.e2 = furman_pivi_surface['e2']
        self.sigmaE = furman_pivi_surface['sigmaE']

        # Parameters for rediffused electrons
        self.p1RInf = furman_pivi_surface['p1RInf']
        self.eR = furman_pivi_surface['eR']
        self.r = furman_pivi_surface['r']
        self.q = furman_pivi_surface['q']
        self.r1 = furman_pivi_surface['r1']
        self.r2 = furman_pivi_surface['r2']

        # Parameters for true secondaries
        self.deltaTSHat = furman_pivi_surface['deltaTSHat']
        self.eHat0 = furman_pivi_surface['eHat0']
        self.s = furman_pivi_surface['s']
        self.t1 = furman_pivi_surface['t1']
        self.t2 = furman_pivi_surface['t2']

        if self.use_ECLOUD_theta0_dependence:
            # Emax(theta) as in ECLOUD module
            self.t3 = 0.7
            self.t4 = 1.
        else:
            self.t3 = furman_pivi_surface['t3']
            self.t4 = furman_pivi_surface['t4']

        if self.exclude_rediffused:
            print(('Secondary emission model: Furman-Pivi excluding rediffused, s=%.4f' % (self.s)))
        else:
            print(('Secondary emission model: Furman-Pivi, s=%.4f' % (self.s)))

    def SEY_model_evol(self, Dt):
        pass

    def SEY_process(self, E_impact_eV, costheta_impact, i_impact):
        """
        Decides event type for each MP colliding with energy E_impact_eV and
        incident angle costheta_impact.
        Returns the SEY components as well as flags defining the event type of
        each MP. Does not rescale the MPs (that is done in impacts_on_surface).
        """
        # Furman-Pivi algorithm
        # (1): Compute emission angles and energy
        # Implemented in the impact_management_class.

        # (2): Compute delta_e, delta_r, delta_ts
        delta_e, delta_r, delta_ts = self.yield_fun_furman_pivi(E_impact_eV, costheta_impact)

        # (3): Generate probability of number of electrons created
        # Implemented in the impact_management_class.

        # Decide on type
        rand = random.rand(E_impact_eV.size)
        if self.exclude_rediffused:
            flag_truesec = rand > delta_e
            flag_backscattered = (~flag_truesec)
            return flag_backscattered, None, flag_truesec, delta_e, None, delta_ts
        else:
            flag_truesec = rand > delta_e + delta_r
            flag_backscattered = np.logical_and(~flag_truesec, rand < delta_e)
            flag_rediffused = np.logical_and(~flag_truesec, ~flag_backscattered)

        # (4): Generate number of secondaries for every impact
        # In impacts_on_surface

        # (5): Delete if n = 0
        # Done automatically by the MP system.

        # (6): Generate energy:
        # In impacts_on_surface

        return flag_backscattered, flag_rediffused, flag_truesec, delta_e, delta_r, delta_ts

    def yield_fun_furman_pivi(self, E, costheta, check=True):
        delta_e = self.delta_e(E, costheta)
        delta_r = self.delta_r(E, costheta)
        delta_ts = self.delta_ts(E, costheta)
        if check and (delta_e + delta_r >= 1).any():
            raise ValueError('delta_e + delta_r is greater than 1')
        return delta_e, delta_r, delta_ts

    def delta_e(self, E_impact_eV, costheta_impact):
        """
        SEY component of backscattered electrons (elastically scattered).
        (25) in FP paper.
        """
        exp_factor = -(np.abs(E_impact_eV - self.eEHat) / self.w)**self.p / self.p
        delta_e0 = self.p1EInf + (self.p1Ehat - self.p1EInf) * np.exp(exp_factor)
        if self.flag_costheta_delta_scale:
            if self.use_ECLOUD_theta0_dependence:
                angular_factor = 1.
            else:
                angular_factor = 1. + self.e1 * (1. - costheta_impact**self.e2)
        else:
            angular_factor = 1

        return delta_e0 * angular_factor

    def delta_r(self, E_impact_eV, costheta_impact):
        """
        SEY component of rediffused electrons (not in ECLOUD model).
        (28) in FP paper.
        """
        exp_factor = -(E_impact_eV / self.eR)**self.r
        delta_r0 = self.p1RInf * (1. - np.exp(exp_factor))
        if self.flag_costheta_delta_scale:
            angular_factor = 1. + self.r1 * (1. - costheta_impact**self.r2)
        else:
            angular_factor = 1

        return delta_r0 * angular_factor

    def delta_ts(self, E_impact_eV, costheta_impact):
        """
        SEY component of true secondaries.
        (31) in FP paper.
        """
        if self.flag_costheta_Emax_shift:
            eHat = self.eHat0 * (1. + self.t3 * (1. - costheta_impact**self.t4))
        else:
            eHat = self.eHat0
        delta_ts0 = self.deltaTSHat * self._D(E_impact_eV / eHat)
        if self.flag_costheta_delta_scale:
            if self.use_ECLOUD_theta0_dependence:
                angular_factor = np.exp(0.5 * (1. - costheta_impact))
            else:
                angular_factor = 1. + self.t1 * (1. - costheta_impact**self.t2)
        else:
            angular_factor = 1

        return delta_ts0 * angular_factor

    def _D(self, x):
        """(32) in FP paper"""
        s = self.s
        return s * x / (s - 1 + x**s)

    def get_energy_backscattered(self, E_0):
        """
        Inverse transform sampling of (26) in the Furman-Pivi paper.
        Returns emission energies for backscattered electrons.
        """
        sqrt2 = np.sqrt(2)
        uu = random.rand(len(E_0))

        if self.use_modified_sigmaE:
            aa = 1.88
            bb = 2.5
            cc = 1e-2
            dd = 1.5e2
            sigmaE_modified = (self.sigmaE - aa) + bb * (1 + np.tanh(cc * (E_0 - dd)))
            return E_0 + sqrt2 * sigmaE_modified * erfinv((uu - 1) * erf(E_0 / (sqrt2 * sigmaE_modified)))
        else:
            return E_0 + sqrt2 * self.sigmaE * erfinv((uu - 1) * erf(E_0 / (sqrt2 * self.sigmaE)))

    def get_energy_rediffused(self, E0):
        """
        Inverse transform sampling of (29) in the Furman-Pivi paper.
        Returns emission energies for rediffused electrons.
        """
        uu = random.rand(len(E0))
        return uu**(1 / (self.q + 1)) * E0

    def _true_sec_energy_CDF(self, nn, energy):
        """
        Gives the value of the CDF corresponding to nn emitted.
        Returns the value of the CDF as well as the area under the PDF before
        normalisation.
        """
        if isinstance(nn, int) or isinstance(nn, np.float64):
            eps_curr = self.eps_n[int(nn - 1)]
            p_n_curr = self.p_n[int(nn - 1)]
        else:
            eps_curr = np.array([self.eps_n[int(ii - 1)] for ii in nn])
            p_n_curr = np.array([self.p_n[int(ii - 1)] for ii in nn])

        cdf = gammainc(p_n_curr, energy / eps_curr)
        area = cdf
        cdf = cdf / area[-1]
        return cdf, area

    def get_energy_true_sec(self, nn, E_0):
        """Returns emission energies for true secondary electrons."""
        if self.use_ECLOUD_energy:
            Ngen = len(E_0)
            return ee.sec_energy_hilleret_model2(self.switch_no_increase_energy, Ngen, self.sigmafit, self.mufit, self.E_th, E_0, self.thresh_low_energy)
        else:
            if len(E_0) == 0:
                return np.array([])
            p_n = self.p_n
            eps_n = self.eps_n

            eps_vec = np.array([eps_n[int(ii - 1)] for ii in nn])
            p_n_vec = np.array([p_n[int(ii - 1)] for ii in nn])
            _, area = self._true_sec_energy_CDF(nn, energy=E_0)  # Putting energy=E_0 gives area under the PDF
            normalisation = 1. / area
            uu = random.rand(len(E_0))
            xx = uu / (normalisation)

            xx[xx < 1e-12] = 0.0  # gammaincinv returns nan if xx is too small but not zero

            return eps_vec * gammaincinv(p_n_vec, xx)

    def inverse_repeat(self, a, repeats, axis):
        """The inverse of numpy.repeat(a, repeats, axis)"""
        if isinstance(repeats, int):
            indices = np.arange(a.shape[axis] / repeats, dtype=np.int) * repeats
        else:  # assume array_like of int
            indices = np.cumsum(repeats) - 1
        return a.take(indices, axis)

    def impacts_on_surface(self, mass, nel_impact, x_impact, y_impact, z_impact,
                           vx_impact, vy_impact, vz_impact, Norm_x, Norm_y, i_found,
                           v_impact_n, E_impact_eV, costheta_impact, nel_mp_th, flag_seg):

        flag_backscattered, flag_rediffused, flag_truesec, \
            delta_e, delta_r, delta_ts = self.SEY_process(E_impact_eV, costheta_impact, i_found)

        nel_replace = nel_impact.copy()
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

        # Backscattered
        if self.use_ECLOUD_energy:
            vx_replace[flag_backscattered], vy_replace[flag_backscattered] = ee.specular_velocity(
                vx_impact[flag_backscattered], vy_impact[flag_backscattered],
                Norm_x[flag_backscattered], Norm_y[flag_backscattered], v_impact_n[flag_backscattered])
        else:
            En_backscattered_eV = self.get_energy_backscattered(E_impact_eV[flag_backscattered])
            N_backscattered = np.sum(flag_backscattered)
            vx_replace[flag_backscattered], vy_replace[flag_backscattered], vz_replace[flag_backscattered] = self.angle_dist_func(
                N_backscattered, En_backscattered_eV, Norm_x[flag_backscattered], Norm_y[flag_backscattered], mass)

        if not self.exclude_rediffused:
            # Rediffused
            En_rediffused_eV = self.get_energy_rediffused(E_impact_eV[flag_rediffused])
            N_rediffused = np.sum(flag_rediffused)
            vx_replace[flag_rediffused], vy_replace[flag_rediffused], vz_replace[flag_rediffused] = self.angle_dist_func(
                N_rediffused, En_rediffused_eV, Norm_x[flag_rediffused], Norm_y[flag_rediffused], mass)

        # True secondary
        N_true_sec = np.sum(flag_truesec)
        n_add_total = 0
        n_emit_truesec_MPs = np.zeros_like(flag_truesec, dtype=int)
        if N_true_sec > 0:
            if self.exclude_rediffused:
                delta_ts_prime = delta_ts[flag_truesec] / (1 - delta_e[flag_truesec])
            else:
                delta_ts_prime = delta_ts[flag_truesec] / (1 - delta_e[flag_truesec] - delta_r[flag_truesec])  # delta_ts^prime in FP paper, eq. (39)

            # Decide how many MPs to be emitted
            n_emit_truesec_MPs[flag_truesec] = random.poisson(lam=delta_ts_prime)  # Using (45)
            n_emit_truesec_MPs_flag_true_sec = n_emit_truesec_MPs[flag_truesec]
            # Cut above M_cut
            flag_above_th = (n_emit_truesec_MPs_flag_true_sec > self.M_cut)
            Nabove_th = np.sum(flag_above_th)
            i_attempt = 0
            while Nabove_th > 0:
                n_emit_truesec_MPs_flag_true_sec[flag_above_th] = random.poisson(delta_ts_prime[flag_above_th])
                if i_attempt>10:
                    n_emit_truesec_MPs_flag_true_sec[flag_above_th] = np.clip(
                            n_emit_truesec_MPs_flag_true_sec[flag_above_th], 0, self.M_cut)
                flag_above_th = (n_emit_truesec_MPs_flag_true_sec > self.M_cut)
                Nabove_th = np.sum(flag_above_th)
                i_attempt += 1
            n_emit_truesec_MPs[flag_truesec] = n_emit_truesec_MPs_flag_true_sec

            # MPs to be replaced
            flag_above_zero = (n_emit_truesec_MPs_flag_true_sec > 0)  # I exclude the absorbed
            flag_truesec_and_above_zero = flag_truesec & (n_emit_truesec_MPs > 0)
            En_truesec_eV = self.get_energy_true_sec(
                nn=n_emit_truesec_MPs_flag_true_sec[flag_above_zero],
                E_0=E_impact_eV[flag_truesec_and_above_zero])  # First generated MPs

            N_true_sec = np.sum(flag_above_zero)

            # Add new MPs
            n_add = n_emit_truesec_MPs - 1
            n_add[n_add < 0] = 0
            n_add_total = np.sum(n_add)
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
                n_emit_truesec_MPs_extended = np.repeat(n_emit_truesec_MPs_flag_true_sec, n_add[flag_truesec])
                En_truesec_eV_add = self.get_energy_true_sec(nn=n_emit_truesec_MPs_extended, E_0=E_impact_eV_add)

                # Ensure energy conservation in each event
                if self.conserve_energy:
                    En_truesec_eV_extended = np.repeat(En_truesec_eV, n_add[flag_truesec][flag_above_zero])  # Energy of replaced MPs
                    E_impact_eV_add = np.repeat(E_impact_eV[flag_truesec][flag_above_zero], n_add[flag_truesec][flag_above_zero])  # Energy of new MPs
                    En_emit_eV_event_add = En_truesec_eV_extended + En_truesec_eV_add  # Total energy emitted in each trusec event

                    flag_violation = (En_emit_eV_event_add > E_impact_eV_add)  # Violation flag for the new MPs
                    flag_violation_replace = self.inverse_repeat(flag_violation, repeats=n_add[flag_truesec][flag_above_zero], axis=None)  # Violation flag for the replaced MPs

                    N_violations = np.sum(flag_violation)
                    while N_violations > 0:
                        # Regenerating energies for the new MPs
                        En_truesec_eV_add[flag_violation] = self.get_energy_true_sec(
                            nn=n_emit_truesec_MPs_extended[flag_violation],
                            E_0=E_impact_eV_add[flag_violation])

                        # Regenerating energies for the replaced MPs
                        En_truesec_eV[flag_violation_replace] = self.get_energy_true_sec(
                            nn=n_emit_truesec_MPs_flag_true_sec[flag_above_zero][flag_violation_replace],
                            E_0=E_impact_eV[flag_truesec][flag_above_zero][flag_violation_replace])

                        En_truesec_eV_extended = np.repeat(En_truesec_eV, n_add[flag_truesec][flag_above_zero])
                        En_emit_eV_event_add = En_truesec_eV_extended + En_truesec_eV_add

                        flag_violation = (En_emit_eV_event_add > E_impact_eV_add)
                        flag_violation_replace = self.inverse_repeat(flag_violation, repeats=n_add[flag_truesec][flag_above_zero], axis=None)

                        N_violations = np.sum(flag_violation)

                # Replace velocities
                vx_replace[flag_truesec_and_above_zero], vy_replace[flag_truesec_and_above_zero], vz_replace[flag_truesec_and_above_zero] = self.angle_dist_func(
                    N_true_sec, En_truesec_eV, Norm_x[flag_truesec_and_above_zero], Norm_y[flag_truesec_and_above_zero], mass)
                # New velocities
                vx_new_MPs, vy_new_MPs, vz_new_MPs = self.angle_dist_func(
                    n_add_total, En_truesec_eV_add, norm_x_add, norm_y_add, mass)

                if flag_seg:
                    i_seg_new_MPs = np.repeat(i_found, n_add)
                else:
                    i_seg_new_MPs = None

            # Handle absorbed MPs
            flag_truesec_and_zero = flag_truesec & (n_emit_truesec_MPs == 0)
            nel_replace[flag_truesec_and_zero] = 0.0
            vx_replace[flag_truesec_and_zero] = 0.0
            vy_replace[flag_truesec_and_zero] = 0.0
            vz_replace[flag_truesec_and_zero] = 0.0
            x_replace[flag_truesec_and_zero] = 0.0
            y_replace[flag_truesec_and_zero] = 0.0
            z_replace[flag_truesec_and_zero] = 0.0

        if n_add_total == 0:
            nel_new_MPs = np.array([])
            x_new_MPs = np.array([])
            y_new_MPs = np.array([])
            z_new_MPs = np.array([])
            vx_new_MPs = np.array([])
            vy_new_MPs = np.array([])
            vz_new_MPs = np.array([])
            i_seg_new_MPs = np.array([])

        # Elastic and rediffused events emit 1 MP
        n_emit_MPs = n_emit_truesec_MPs
        n_emit_MPs[flag_backscattered] = 1
        if not self.exclude_rediffused:
            n_emit_MPs[flag_rediffused] = 1

        nel_emit_tot_events = nel_impact * n_emit_MPs
        events = flag_truesec.astype(int)
        events[n_emit_MPs == 0] = 3  # Absorbed MPs

        if self.exclude_rediffused:
            pass
        else:
            events = events + 2 * flag_rediffused.astype(int)
        event_type = events

        # extended_event_type keeps track of the event type for new MPs, it is
        # needed for the extraction of emission-energy distributions.
        if n_add_total != 0:
            events_add = np.repeat(events, n_add)
            events = np.concatenate([events, events_add])
        extended_event_type = events

        event_info = {'extended_event_type': extended_event_type,
                      }

        return nel_emit_tot_events, event_type, event_info,\
            nel_replace, x_replace, y_replace, z_replace, vx_replace, vy_replace, vz_replace, i_seg_replace,\
            nel_new_MPs, x_new_MPs, y_new_MPs, z_new_MPs, vx_new_MPs, vy_new_MPs, vz_new_MPs, i_seg_new_MPs

    ############################################################################
    #  The following functions are not used in the simulation code but are     #
    #  provided here for use in tests and development.                         #
    ############################################################################
    def backscattered_energy_PDF(self, energy, E_0):
        """The PDF for backscattered electrons."""
        if self.use_modified_sigmaE:
            aa = 1.88
            bb = 2.5
            cc = 1e-2
            dd = 1.5e2
            sigmaE_modified = (self.sigmaE - aa) + bb * (1 + np.tanh(cc * (E_0 - dd)))
            sigma_e = sigmaE_modified
        else:
            sigma_e = self.sigmaE
        ene = energy - E_0
        a = 2 * np.exp(-(ene)**2 / (2 * sigma_e**2))
        c = (np.sqrt(2 * np.pi) * sigma_e * erf(E_0 / (np.sqrt(2) * sigma_e)))
        return a / c

    def backscattered_energy_CDF(self, energy, E_0, sigma_e=2.):
        """The CDF for backscattered electrons."""
        sqrt2 = np.sqrt(2)
        return 1 - erf((E_0 - energy) / (sqrt2 * sigma_e)) / erf(E_0 / (sqrt2 * sigma_e))

    def rediffused_energy_PDF(self, energy, E_0, qq=0.5):
        for ene in E_0:
            if ene < 0:
                raise ValueError('Impacting energy E_0 cannot be negative')
        prob_density = (qq + 1) * energy**qq / E_0**(qq + 1)
        return prob_density

    def rediffused_energy_CDF(self, energy, E_0, qq=0.5):
        return energy**(qq + 1) / E_0**(qq + 1)

    def true_sec_energy_PDF(self, delta_ts, nn, E_0, energy):
        """The PDF for true secondary electrons."""
        if nn == 0:
            raise ValueError('nn = 0, you cannot emit zero electrons.')
        p_n = self.p_n
        eps_n = self.eps_n

        nn_all = np.arange(0, self.M_cut + 1, 1)

        if self.choice == 'poisson':
            P_n_ts = np.squeeze(delta_ts**nn_all / factorial(nn_all) * np.exp(-delta_ts))
        elif self.choice == 'binomial':
            p = delta_ts / self.M_cut
            P_n_ts = np.squeeze(binom(self.M_cut, nn) * (p)**nn_all * (1 - p)**(self.M_cut - nn_all))
        else:
            raise ValueError('choice must be either \'poisson\' or \'binomial\'')

        P_n_ts = P_n_ts / np.sum(P_n_ts)
        P_n_ts_return = P_n_ts[int(nn)]
        eps_curr = eps_n[int(nn - 1)]
        p_n_curr = p_n[int(nn - 1)]

        if E_0 == 0:
            F_n = 0
        else:
            F_n = 1
        f_n_ts = F_n * energy**(p_n_curr - 1) * np.exp(-energy / eps_curr)
        area = scipy.integrate.simps(f_n_ts, energy)
        f_n_ts = f_n_ts / area  # normalisation

        return f_n_ts, P_n_ts_return

    def average_true_sec_energy_PDF(self, delta_ts, E_0, energy):
        nns = np.arange(1, self.M_cut + 1, 1)
        average_f_n_ts = np.zeros_like(energy)
        for ii in nns:
            f_n_ts, P_n_ts = self.true_sec_energy_PDF(delta_ts=delta_ts, nn=ii, E_0=E_0, energy=energy)
            if False:#self.conserve_energy:
                p_n_curr = self.p_n[ii - 1]
                eps_n_curr = self.eps_n[ii - 1]
                factor = eps_n_curr**p_n_curr * gamma(p_n_curr) * gammainc(ii * p_n_curr, E_0 / eps_n_curr)
                factor = 1.
                if ii - 1 == 0:
                    average_f_n_ts = average_f_n_ts + f_n_ts * P_n_ts * ii * 1. / factor
                else:
                    average_f_n_ts = average_f_n_ts + f_n_ts * P_n_ts * ii * gammainc((ii - 1) * p_n_curr, (E_0 - energy) / eps_n_curr) / factor
            else:
                average_f_n_ts = average_f_n_ts + f_n_ts * P_n_ts * ii
        area = scipy.integrate.simps(average_f_n_ts, energy)
        return average_f_n_ts / area
    ############################################################################
    ############################################################################

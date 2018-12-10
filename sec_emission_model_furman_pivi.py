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
import scipy
from scipy.special import gamma
from scipy.special import gammainc
from scipy.special import gammaincinv
from scipy.special import binom
from scipy.special import erf
from scipy.special import erfinv
from scipy.misc import factorial
from scipy.integrate import cumtrapz

# def inverse_CDF_ts_energy(delta_e, delta_r, delta_ts, n, pn, epsn):
#     delta_prime_ts = delta_ts / (1. - delta_e - delta_r)
#     PnTS_prime = delta_prime_ts**n / np.math.factorial(n) * np.exp(-delta_prime_ts)
#     Fn = PnTS_prime / ()(eps[n-1]**pn[n-1] * gamma(pn[n-1]))**n * )


class SEY_model_furman_pivi():

    event_types = {
        0: 'elast',
        1: 'true',
        2: 'rediff',
        3: 'absorb',
    }

    def __init__(self,  # eHat0, deltaTSHat,
                 E_th=None, sigmafit=None, mufit=None,
                 switch_no_increase_energy=0, thresh_low_energy=None, secondary_angle_distribution=None,
                 s=1.54,
                 M=10
                 ):

        self.M = M
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

        print('Secondary emission model: Furman-Pivi s=%.4f' % (self.s))

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

        flag_truesec = rand > delta_e + delta_r
        flag_backscattered = np.logical_and(~flag_truesec, rand < delta_e)
        flag_rediffused = np.logical_and(~flag_truesec, ~flag_backscattered)

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

        nel_emit = nel_impact # * delta

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

    def backscattered_energy_PDF(self, energy, E_0, sigma_e=2.):
        ene = energy - E_0
        a = 2 * np.exp(-(ene)**2 / (2 * sigma_e**2))
        c = (np.sqrt(2 * np.pi) * sigma_e * erf(E_0 / (np.sqrt(2) * sigma_e)))
        return a / c

    def backscattered_energy_CDF(self, energy, E_0, sigma_e=2.):
        sqrt2 = np.sqrt(2)
        return 1 - erf((E_0 - energy) / (sqrt2 * sigma_e)) / erf(E_0 / (sqrt2 * sigma_e))

    def get_energy_backscattered(self, E_0):
        sqrt2 = np.sqrt(2)
        uu = random.rand(len(E_0))
        return E_0 - sqrt2 * self.sigmaE * erfinv(-(uu - 1) * erf(E_0 / (sqrt2 * self.sigmaE)))  # Inverse transform sampling of (26) in FP paper

    def rediffused_energy_PDF(self, energy, E_0, qq=0.5):
        for ene in E_0:
            if ene < 0:
                raise ValueError('Impacting energy E_0 cannot be negative')
        prob_density = (qq + 1) * energy**qq / E_0**(qq + 1)
        return prob_density

    def rediffused_energy_CDF(self, energy, E_0, qq=0.5):
        return energy**(qq + 1) / E_0**(qq + 1)

    def get_energy_rediffused(self, E0):
        uu = random.rand(len(E0))
        return uu**(1 / (self.q + 1)) * E0  # Inverse transform sampling of (29) in FP paper

    def true_sec_energy_PDF(self, delta_ts, nn, E_0, energy=np.linspace(0.001, 300, num=int(1e5)), choice='poisson', M=10):
        """
        A simplified version of the energy distribution for secondary electrons in
        the Furman-Pivi model. The 'choice' parameter decides wheter to use a poisson
        or a Poisson distribution for the probabilities P_n_ts.
        """
        p_n = self.p_n
        eps_n = self.eps_n

        if choice == 'poisson':
            P_n_ts = np.squeeze(delta_ts / factorial(nn) * np.exp(-delta_ts))
        elif choice == 'binomial':
            p = delta_ts / M
            P_n_ts = np.squeeze(binom(M, nn) * (p)**nn * (1 - p)**(M - nn))
        else:
            raise ValueError('choice must be either \'poisson\' or \'binomial\'')
        P_n_ts_return = P_n_ts

        eps_curr = eps_n[int(nn - 1)]
        p_n_curr = p_n[int(nn - 1)]

        # eps_curr = []
        # p_n_curr = []
        # if type(nn) == int or nn.shape == np.empty(None).shape:
        #     eps_curr = eps_n[int(nn - 1)]
        #     p_n_curr = p_n[int(nn - 1)]
        # else:
        #     for kk in nn:
        #         eps_curr.append(eps_n[int(kk - 1)])
        #         p_n_curr.append(p_n[int(kk - 1)])
        #     eps_curr = np.array(eps_curr)
        #     p_n_curr = np.array(p_n_curr)
        P_n_ts = 1.
        if E_0 == 0:
            F_n = 0
        else:
            F_n = np.squeeze((P_n_ts / ((eps_curr**p_n_curr * gamma(p_n_curr))**nn * gammainc(nn * p_n_curr, E_0 / eps_curr)))**(1. / nn))
        f_n_ts = F_n * energy**(p_n_curr - 1) * np.exp(-energy / eps_curr)
        # area = scipy.integrate.simps(f_n_ts, energy)
        # if area != 0:
        #     f_n_ts = f_n_ts / area  # Normalization

        return f_n_ts, P_n_ts_return

    def true_sec_energy_CDF(self, delta_ts, nn, E_0, energy=np.linspace(0.001, 300, num=int(1e5)), choice='poisson', M=10):
        p_n = self.p_n
        eps_n = self.eps_n

        P_n_ts = 1.
        eps_curr = eps_n[int(nn - 1)]
        p_n_curr = p_n[int(nn - 1)]

        if E_0 == 0:
            F_n = 0
        else:
            F_n = np.squeeze((P_n_ts / ((eps_curr**p_n_curr * gamma(p_n_curr))**nn * gammainc(nn * p_n_curr, E_0 / eps_curr)))**(1. / nn))
        return eps_curr**p_n_curr * F_n * gamma(p_n_curr) * gammainc(p_n_curr, energy / eps_curr)

    # def true_sec_energy_CDF(self, delta_ts, nn, E_0, energy=np.linspace(0.001, 300, num=int(1e5)), choice='poisson', M=10):
    #     f_n_ts, _ = self.true_sec_energy_PDF(delta_ts=delta_ts, nn=nn, E_0=E_0, energy=energy, choice=choice, M=M)
    #     CDF = cumtrapz(f_n_ts, energy, initial=0)
    #     return CDF

    def get_energy_true_sec(self, delta_ts, nn, E_0, choice='poisson', M=10):
        p_n = self.p_n
        eps_n = self.eps_n

        P_n_ts = 1.
        eps_vec = np.array([eps_n[int(ii - 1)] for ii in nn])
        p_n_vec = np.array([p_n[int(ii - 1)] for ii in nn])
        F_n_vec = (P_n_ts / ((eps_vec**p_n_vec * gamma(p_n_vec))**nn * gammainc(nn * p_n_vec, E_0 / eps_vec)))**(1. / nn)
        uu = random.rand(len(E_0))
        xx = uu / (F_n_vec * eps_vec**p_n_vec * gamma(p_n_vec))
        xx[xx < 1e-12] = 0.0  # gammaincinv returns nan if xx is too small but not zero

        return eps_vec * gammaincinv(p_n_vec, xx)

    # def get_energy_true_sec(self, delta_ts, nn, E_0, choice='poisson', M=10):
    #     energy = np.linspace(1e-10, 300, num=int(len(E_0)))
    #     uu = random.rand(len(energy))
    #     out_array = np.empty(1)
    #     if type(nn) is int:
    #         CDF = self.true_sec_energy_CDF(delta_ts=delta_ts, nn=nn, E_0=E_0, choice=choice, energy=energy, M=M)
    #         return np.interp(uu, CDF, energy)
    #     else:
    #         for ii, kk in enumerate(nn):
    #             CDF = self.true_sec_energy_CDF(delta_ts=delta_ts[ii], nn=kk, E_0=E_0[ii], choice=choice, energy=energy, M=M)
    #             out_array = np.concatenate([out_array, np.array([np.interp(uu[ii], CDF, energy)])])
    #         out_array = np.delete(out_array, 0)
    #         return out_array

    def average_true_sec_energy_PDF(self, delta_ts, E_0, energy=np.linspace(0.001, 300, num=int(1e5)), choice='poisson'):
        nns = np.arange(1, 11, 1)
        average_f_n_ts = np.zeros_like(energy)
        for ii in nns:
            f_n_ts, P_n_ts = self.true_sec_energy_PDF(delta_ts=delta_ts, nn=ii, E_0=E_0, choice=choice, energy=energy)
            average_f_n_ts = average_f_n_ts + f_n_ts * P_n_ts
        area = scipy.integrate.simps(average_f_n_ts, energy)
        normalization_constant = 1. / area
        return normalization_constant * average_f_n_ts

    def average_true_sec_energy_CDF(self, delta_ts, E_0, choice='poisson', energy=np.linspace(0.001, 300, num=int(1e5))):
        pdf = self.average_true_sec_energy_PDF(delta_ts=delta_ts, E_0=E_0, choice=choice, energy=energy)
        CDF = cumtrapz(pdf, energy, initial=0)
        return CDF

    # numpy.interp(x, xp, fp, left=None, right=None, period=None)
    def get_energy_average_true_sec(self, delta_ts, E_0, choice='poisson', energy=np.linspace(0.001, 300, num=int(1e5))):
        uu = random.rand(len(delta_ts))
        out_array = np.empty(1)
        if type(delta_ts) is int:
            CDF = self.average_true_sec_energy_CDF(delta_ts=delta_ts, E_0=E_0, choice=choice, energy=energy)
            return np.interp(uu, CDF, energy)
        else:
            for ii, kk in enumerate(delta_ts):
                CDF = self.average_true_sec_energy_CDF(delta_ts=delta_ts[ii], E_0=E_0[ii], choice=choice, energy=energy)
                out_array = np.concatenate([out_array, np.array([np.interp(uu[ii], CDF, energy)])])
        return np.interp(uu, CDF, energy)

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

        # Backscattered
        En_backscattered_eV = self.get_energy_backscattered(E_impact_eV[flag_backscattered])
        N_rediffused = np.sum(flag_backscattered)
        vx_replace[flag_backscattered], vy_replace[flag_backscattered], vz_replace[flag_backscattered] = self.angle_dist_func(
            N_rediffused, En_backscattered_eV, Norm_x[flag_backscattered], Norm_y[flag_backscattered], mass)

        # Rediffused
        En_rediffused_eV = self.get_energy_rediffused(E_impact_eV[flag_rediffused])
        N_rediffused = np.sum(flag_rediffused)
        vx_replace[flag_rediffused], vy_replace[flag_rediffused], vz_replace[flag_rediffused] = self.angle_dist_func(
            N_rediffused, En_rediffused_eV, Norm_x[flag_rediffused], Norm_y[flag_rediffused], mass)

        # True secondary
        N_true_sec = np.sum(flag_truesec)
        n_add_total = 0
        if N_true_sec > 0:
            # delta_ts = self._delta_ts(E_impact_eV[flag_truesec], costheta_impact[flag_truesec])
            delta_e, delta_r, delta_ts = self._yield_fun_furman_pivi(E_impact_eV[flag_truesec], costheta_impact[flag_truesec])
            delta_ts = delta_ts / (1 - delta_e - delta_r)  # delta_ts^prime in FP paper, eq. (39)
            n_add = np.zeros_like(flag_truesec, dtype=int)
            n_add[flag_truesec] = random.poisson(lam=delta_ts)
            n_add_flag_true_sec = n_add[flag_truesec]
            # Cut above M
            flag_above_th = (n_add[flag_truesec] > self.M)
            Nabove_th = np.sum(flag_above_th)
            while Nabove_th > 0:
                n_add_flag_true_sec[flag_above_th] = random.poisson(delta_ts[flag_above_th])
                flag_above_th = (n_add_flag_true_sec > self.M)
                Nabove_th = np.sum(flag_above_th)
            n_add[flag_truesec] = n_add_flag_true_sec

            # MPs to be replaced
            flag_above_zero = (n_add_flag_true_sec > 0)
            flag_truesec_and_above_zero = flag_truesec & (n_add > 0)
            En_truesec_eV = self.get_energy_true_sec(delta_ts=delta_ts[flag_above_zero], nn=n_add_flag_true_sec[flag_above_zero], E_0=E_impact_eV[flag_truesec_and_above_zero], M=self.M)

            N_true_sec = np.sum(flag_above_zero)
            vx_replace[flag_truesec_and_above_zero], vy_replace[flag_truesec_and_above_zero], vz_replace[flag_truesec_and_above_zero] = self.angle_dist_func(
                N_true_sec, En_truesec_eV, Norm_x[flag_truesec_and_above_zero], Norm_y[flag_truesec_and_above_zero], mass)

            flag_truesec_and_zero = flag_truesec & (n_add == 0)
            nel_replace[flag_truesec_and_zero] = 0.0

            # Add new MPs
            clone_idxs = n_add - 1
            clone_idxs[clone_idxs < 0] = 0
            n_add_total = np.sum(clone_idxs)
            if n_add_total != 0:
                # Clone MPs
                x_new_MPs = np.repeat(x_impact, clone_idxs)
                y_new_MPs = np.repeat(y_impact, clone_idxs)
                z_new_MPs = np.repeat(z_impact, clone_idxs)
                norm_x_add = np.repeat(Norm_x, clone_idxs)
                norm_y_add = np.repeat(Norm_y, clone_idxs)
                nel_new_MPs = np.repeat(nel_replace, clone_idxs)
                E_impact_eV_add = np.repeat(E_impact_eV, clone_idxs)

                # Generate new MP properties, angles and energies
                flag_above_one = (n_add[flag_truesec] > 1)
                clone_above_one = clone_idxs[flag_truesec][flag_above_one]
                n_add_extended = np.repeat(n_add_flag_true_sec[flag_above_one], clone_above_one)
                delta_ts_extended = np.repeat(delta_ts[flag_above_one], clone_above_one)

                En_truesec_eV_add = self.get_energy_true_sec(delta_ts=delta_ts_extended, nn=n_add_extended, E_0=E_impact_eV_add, M=self.M)

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

        # Add new MPs to nel_emit_tot_events
        nel_emit_tot_events = np.concatenate([nel_emit_tot_events, nel_new_MPs])

        events = flag_truesec.astype(int)
        if N_true_sec > 0:
            events[n_add == 0] = 3  # Absorbed MPs

        events = events - flag_rediffused.astype(int) - flag_backscattered.astype(int) * 3
        if n_add_total != 0:
            events_add = np.repeat(events, clone_idxs)
            events = np.concatenate([events, events_add])
        event_type = events

        event_info = {}

        return nel_emit_tot_events, event_type, event_info,\
            nel_replace, x_replace, y_replace, z_replace, vx_replace, vy_replace, vz_replace, i_seg_replace,\
            nel_new_MPs, x_new_MPs, y_new_MPs, z_new_MPs, vx_new_MPs, vy_new_MPs, vz_new_MPs, i_seg_new_MPs


class SEY_model_FP_Cu(SEY_model_furman_pivi):

    p_n = np.array([2.5, 3.3, 2.5, 2.5, 2.8, 1.3, 1.5, 1.5, 1.5, 1.5])
    eps_n = np.array([1.5, 1.75, 1, 3.75, 8.5, 11.5, 2.5, 3, 2.5, 3])

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

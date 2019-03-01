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
from sec_emission_model_furman_pivi import SEY_model_furman_pivi


class SEY_model_furman_pivi_variable_MP(SEY_model_furman_pivi):

    event_types = {
        0: 'elast',
        1: 'true',
        2: 'rediff',
    }

    def __init__(self, furman_pivi_surface,
                 E_th=None, sigmafit=None, mufit=None,
                 switch_no_increase_energy=0, thresh_low_energy=None, secondary_angle_distribution=None,
                 ):

        SEY_model_furman_pivi.__init__(self, furman_pivi_surface,
                                       E_th=E_th, sigmafit=sigmafit, mufit=mufit,
                                       switch_no_increase_energy=switch_no_increase_energy,
                                       thresh_low_energy=thresh_low_energy,
                                       secondary_angle_distribution=secondary_angle_distribution,
                                       )
        self.variable_MP_size = True

        print('Secondary emission model: Furman-Pivi, varaible MP size, s=%.4f' % (self.s))

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
        nel_emit = nel_impact * delta

        return nel_emit, flag_backscattered, flag_rediffused, flag_truesec, delta_e, delta_r, delta_ts

    def impacts_on_surface(self, mass, nel_impact, x_impact, y_impact, z_impact,
                           vx_impact, vy_impact, vz_impact, Norm_x, Norm_y, i_found,
                           v_impact_n, E_impact_eV, costheta_impact, nel_mp_th, flag_seg):

        nel_emit_tot_events, flag_backscattered, flag_rediffused, flag_truesec, delta_e, delta_r, delta_ts = self.SEY_process(nel_impact, E_impact_eV, costheta_impact, i_found)

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
            n_add = np.zeros_like(flag_truesec, dtype=int)
            n_add[flag_truesec] = np.ceil(nel_replace[flag_truesec] / nel_mp_th) - 1
            n_add[n_add < 0] = 0.  # in case of underflow
            nel_replace[flag_truesec] = nel_replace[flag_truesec] / (n_add[flag_truesec] + 1.)

            n_add_total = np.sum(n_add)

            # MPs to be replaced
            # En_truesec_eV = ee.sec_energy_hilleret_model2(
            #     self.switch_no_increase_energy, N_true_sec, self.sigmafit, self.mufit,
            #     self.E_th, E_impact_eV[flag_truesec], self.thresh_low_energy)
            delta_ts_prime = delta_ts / (1. - delta_e - delta_r)
            En_truesec_eV = self.get_energy_average_true_sec(delta_ts=delta_ts_prime[flag_truesec], E_0=E_impact_eV[flag_truesec])

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
                delta_ts_prime_add = np.repeat(delta_ts_prime, n_add)

                # Generate new MP properties, angles and energies
                # En_truesec_eV_add = ee.sec_energy_hilleret_model2(
                #     self.switch_no_increase_energy, n_add_total, self.sigmafit, self.mufit,
                #     self.E_th, E_impact_eV_add, self.thresh_low_energy)
                En_truesec_eV_add = self.get_energy_average_true_sec(delta_ts=delta_ts_prime_add, E_0=E_impact_eV_add)

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

        # extended_nel_emit_tot_events used for extraction of energy distributions
        extended_nel_emit_tot_events = np.concatenate([nel_replace, nel_new_MPs])

        events = flag_truesec.astype(int)

        events = events + 2 * flag_rediffused.astype(int)
        event_type = events
        if n_add_total != 0:
            events_add = np.repeat(events, n_add)
            events = np.concatenate([events, events_add])
        extended_event_type = events

        event_info = {'extended_nel_emit_tot_events': extended_nel_emit_tot_events,
                      'extended_event_type': extended_event_type,
                      }

        return nel_emit_tot_events, event_type, event_info,\
            nel_replace, x_replace, y_replace, z_replace, vx_replace, vy_replace, vz_replace, i_seg_replace,\
            nel_new_MPs, x_new_MPs, y_new_MPs, z_new_MPs, vx_new_MPs, vy_new_MPs, vz_new_MPs, i_seg_new_MPs

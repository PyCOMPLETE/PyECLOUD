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

import hist_for as histf
import numpy as np
from scipy.constants import e as qe
from impact_management_class import impact_management


class impact_management_perfect_absorber(impact_management):

    #@profile
    def backtrack_and_second_emiss(self, old_pos, MP_e):

        self.Nel_impact_last_step = 0.
        self.Nel_emit_last_step = 0.
        self.En_imp_last_step_eV = 0.
        self.En_emit_last_step_eV = 0.

        if MP_e.N_mp > 0:

            x_mp_old = old_pos.x_mp
            y_mp_old = old_pos.y_mp
            z_mp_old = old_pos.z_mp
            x_mp = MP_e.x_mp
            y_mp = MP_e.y_mp
            z_mp = MP_e.z_mp
            vx_mp = MP_e.vx_mp
            vy_mp = MP_e.vy_mp
            vz_mp = MP_e.vz_mp
            nel_mp = MP_e.nel_mp
            N_mp = MP_e.N_mp
            chamb = self.chamb
            bias_x_hist = self.bias_x_hist
            Dx_hist = self.Dx_hist
            En_hist_max = self.En_hist_max
            DEn_hist = self.DEn_hist
            flag_seg = self.flag_seg
            scrub_en_th = self.scrub_en_th

            ## impact management
            N_mp_old = N_mp

            #~ # check flag_impact array has right size (if not regenerate it)
            #~ if len(self.flag_impact) != len(x_mp):
            #~     self.flag_impact = np.array(len(x_mp)*[False])

            # reset flag_impact array
            self.flag_impact = np.zeros_like(x_mp, dtype=bool)

            # detect impact
            self.flag_impact[:N_mp] = chamb.is_outside(x_mp[0:N_mp], y_mp[0:N_mp])  # (((x_mp[0:N_mp]/x_aper)**2 + (y_mp[0:N_mp]/y_aper)**2)>=1);

            Nimpact = int(np.sum(self.flag_impact))

            if Nimpact > 0:

                if flag_seg:
                    raise ValueError('Segment identification not implemented for perfect absorber! Sorry...')

                # load segment endpoints
                x_in = x_mp_old[self.flag_impact[:N_mp_old]]; y_in = y_mp_old[self.flag_impact[:N_mp_old]]; z_in = z_mp_old[self.flag_impact[:N_mp_old]]

                # just to have them back in the chamber
                x_emit = x_in
                y_emit = y_in
                z_emit = z_in

                # load velocities and charges
                vx_impact = vx_mp[self.flag_impact]; vy_impact = vy_mp[self.flag_impact]; vz_impact = vz_mp[self.flag_impact]
                nel_impact = nel_mp[self.flag_impact]

                # compute impact velocities, energy and angle
                v_impact_mod = np.sqrt(vx_impact * vx_impact + vy_impact * vy_impact + vz_impact * vz_impact)
                E_impact_eV = 0.5 * MP_e.mass / qe * v_impact_mod * v_impact_mod

                #electron histogram
                histf.compute_hist(x_emit, nel_impact, bias_x_hist, Dx_hist, self.nel_impact_hist_tot)
                histf.compute_hist(x_emit, nel_impact * (E_impact_eV > scrub_en_th), bias_x_hist, Dx_hist, self.nel_impact_hist_scrub)
                histf.compute_hist(x_emit, nel_impact * E_impact_eV, bias_x_hist, Dx_hist, self.energ_eV_impact_hist)

                En_imp_hist = E_impact_eV.copy()
                En_imp_hist[En_imp_hist > En_hist_max] = En_hist_max
                histf.compute_hist(En_imp_hist, nel_impact, 0., DEn_hist, self.En_hist_line)

                nel_emit = 0 * nel_impact

                self.Nel_impact_last_step = np.sum(nel_impact)
                self.Nel_emit_last_step = np.sum(nel_emit)
                self.En_imp_last_step_eV = np.sum(E_impact_eV * nel_impact)

                # absorb
                x_mp[self.flag_impact] = x_emit
                y_mp[self.flag_impact] = y_emit
                z_mp[self.flag_impact] = z_emit
                vx_mp[self.flag_impact] = 0.
                vy_mp[self.flag_impact] = 0.
                vz_mp[self.flag_impact] = 0.
                nel_mp[self.flag_impact] = 0.

                MP_e.x_mp = x_mp
                MP_e.y_mp = y_mp
                MP_e.z_mp = z_mp
                MP_e.vx_mp = vx_mp
                MP_e.vy_mp = vy_mp
                MP_e.vz_mp = vz_mp
                MP_e.N_mp = N_mp

        return MP_e

    def extract_sey_curves(self, n_rep, E_impact_eV_test, cos_theta_test, charge, mass):

        del_true_mat = np.zeros((len(cos_theta_test), len(E_impact_eV_test)))
        del_elast_mat = np.zeros((len(cos_theta_test), len(E_impact_eV_test)))
        print('Extracting SEY curves...')
        print(':-P')
        print('Done extracting SEY curves.')

        return del_true_mat, del_elast_mat

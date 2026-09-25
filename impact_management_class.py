#-Begin-preamble-------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 8.7.1
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
from . import hist_for as histf
from . import seg_impact as segi
from scipy.constants import e as qe

class impact_management(object):
    def __init__(
            self, chamb, sey_mod,
            Dx_hist, scrub_en_th, Nbin_En_hist,
            En_hist_max, Nbin_lifetime_hist=None,
            lifetime_hist_max=None, flag_lifetime_hist=False,
            flag_seg=False, flag_En_hist_seg=False,
            cos_angle_width=0.05, flag_cos_angle_hist=True):

        print('Start impact man. init.')

        if flag_seg and not(hasattr(chamb, 'N_vert')):
            raise ValueError(
                """flag_seg can be True only with polygonal chambers!!!!""")

        self.chamb = chamb
        self.sey_mod = sey_mod
        self.Dx_hist = Dx_hist
        self.scrub_en_th = scrub_en_th
        self.Nbin_En_hist = Nbin_En_hist
        self.En_hist_max = En_hist_max
        self.flag_seg = flag_seg
        self.flag_En_hist_seg = flag_En_hist_seg

        xg_hist = np.arange(0, chamb.x_aper + 2. * Dx_hist, Dx_hist, float)
        xgr_hist = xg_hist[1:]
        xgr_hist = xgr_hist[::-1]  # reverse array
        xg_hist = np.concatenate((-xgr_hist, xg_hist), 0)
        Nxg_hist = len(xg_hist)
        bias_x_hist = np.min(xg_hist)

        self.En_g_hist = np.linspace(
            0., En_hist_max, Nbin_En_hist)  # hist. grid
        self.DEn_hist = self.En_g_hist[1] - self.En_g_hist[0]  # hist. step

        self.flag_cos_angle_hist = flag_cos_angle_hist
        if flag_cos_angle_hist:
            self.cos_angle_width = cos_angle_width
            N_angles = int(1. / cos_angle_width) + 1
            self.cos_angle_hist = np.zeros(N_angles, float)
            print('Saving cosine of angle of incident electrons.')
        else:
            print('Not saving cosine of angle of incident electrons.')

        self.xg_hist = xg_hist
        self.Nxg_hist = Nxg_hist
        self.bias_x_hist = bias_x_hist

        self.Nel_impact_last_step = None
        self.Nel_emit_last_step = None
        self.En_imp_last_step_eV = None
        self.En_emit_last_step_eV = None

        self.nel_impact_hist_tot = np.zeros(Nxg_hist, float)
        self.nel_impact_hist_scrub = np.zeros(Nxg_hist, float)
        self.energ_eV_impact_hist = np.zeros(Nxg_hist, float)
        self.En_hist_line = np.zeros(Nbin_En_hist, float)

        self.flag_lifetime_hist = flag_lifetime_hist

        if flag_lifetime_hist:
            self.Nbin_lifetime_hist = Nbin_lifetime_hist
            self.lifetime_hist_max = lifetime_hist_max
            self.lifetime_g_hist = np.linspace(
                0., lifetime_hist_max, Nbin_lifetime_hist)  # hist. grid
            # hist. step
            self.Dt_lifetime_hist = self.lifetime_g_hist[1] - \
                self.lifetime_g_hist[0]
            self.lifetime_hist_line = np.zeros(Nbin_lifetime_hist, float)

        if flag_seg:
            self.nel_hist_impact_seg = np.zeros(chamb.N_vert, float)
            self.nel_hist_emit_seg = np.zeros(chamb.N_vert, float)
            self.energ_eV_impact_seg = np.zeros(chamb.N_vert, float)
            if flag_En_hist_seg:
                self.seg_En_hist_lines = [
                    np.zeros(Nbin_En_hist, float) for _ in range(chamb.N_vert)]

        print('Done impact man. init.')

    def reset_impact_hist_tot(self):
        self.nel_impact_hist_tot *= 0.

    def reset_impact_hist_scrub(self):
        self.nel_impact_hist_scrub *= 0.

    def reset_energ_eV_impact_hist(self):
        self.energ_eV_impact_hist *= 0.

    def reset_En_hist_line(self):
        self.En_hist_line *= 0.

    def reset_seg_En_hist_lines(self):
        for ii in range(self.chamb.N_vert):
            self.seg_En_hist_lines[ii] *= 0.

    def reset_hist_impact_seg(self):
        if self.flag_seg:
            self.nel_hist_impact_seg *= 0.

    def reset_hist_emit_seg(self):
        if self.flag_seg:
            self.nel_hist_emit_seg *= 0.

    def reset_energ_impact_seg(self):
        if self.flag_seg:
            self.energ_eV_impact_seg *= 0.

    def reset_cos_angle_hist(self):
        self.cos_angle_hist *= 0

    def reset_lifetime_hist_line(self):
        self.lifetime_hist_line *= 0.

    def backtrack_and_second_emiss(self, old_pos, MP_e, tt_curr=None):

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

            N_mp_old = MP_e.N_mp

            nel_mp_th = MP_e.nel_mp_split
            chamb = self.chamb
            sey_mod = self.sey_mod
            bias_x_hist = self.bias_x_hist
            Dx_hist = self.Dx_hist
            En_hist_max = self.En_hist_max
            DEn_hist = self.DEn_hist

            if self.flag_lifetime_hist:
                Dt_lifetime_hist = self.Dt_lifetime_hist

            flag_seg = self.flag_seg
            scrub_en_th = self.scrub_en_th

            # impact management

            flag_impact = np.zeros_like(x_mp, dtype=bool)
            self.flag_impact = flag_impact

            # detect impact
            flag_impact[:N_mp_old] = chamb.is_outside(
                x_mp[0:N_mp_old], y_mp[0:N_mp_old])
            
            idx = np.flatnonzero(flag_impact)
            Nimpact = int(np.sum(flag_impact))

            if Nimpact > 0:

                # load segment endpoints
                x_in = x_mp_old[idx]
                y_in = y_mp_old[idx]
                z_in = z_mp_old[idx]

                x_out = x_mp[idx]
                y_out = y_mp[idx]
                z_out = z_mp[idx]

                # backtracking and surface normal generation
                [x_impact, y_impact, z_impact, Norm_x, Norm_y, i_found] =\
                    chamb.impact_point_and_normal(
                        x_in, y_in, z_in, x_out, y_out, z_out)

                # load velocities and charges
                vx_impact = vx_mp[idx]
                vy_impact = vy_mp[idx]
                vz_impact = vz_mp[idx]
                nel_impact = nel_mp[idx]

                # add to lifetime histogram
                if self.flag_lifetime_hist:
                    lifetime_impact = tt_curr - MP_e.t_last_impact[idx]
                    if sum(MP_e.t_last_impact[idx] > 0) > 0:
                        histf.compute_hist(lifetime_impact[MP_e.t_last_impact[idx] > 0],
                                           nel_impact[MP_e.t_last_impact[idx] > 0], 
                                           0., Dt_lifetime_hist, self.lifetime_hist_line)

                    MP_e.t_last_impact[idx] = tt_curr

                # compute impact velocities, energy and angle
                v_impact_mod = np.sqrt(
                    vx_impact * vx_impact + vy_impact * vy_impact + vz_impact * vz_impact)
                E_impact_eV = 0.5 * MP_e.mass / qe * v_impact_mod * v_impact_mod
                v_impact_n = vx_impact * Norm_x + vy_impact * Norm_y
                # Use np.abs to rule out negative values, which can happen in very seldom fringe cases.
                # Mathematically correct would be -(v_impact_n)/v_impact_mod
                costheta_impact = np.abs(v_impact_n / v_impact_mod)

                # electron histogram
                histf.compute_hist(
                    x_impact, nel_impact, bias_x_hist, Dx_hist, self.nel_impact_hist_tot)
                histf.compute_hist(x_impact, nel_impact * (E_impact_eV > scrub_en_th),
                                   bias_x_hist, Dx_hist, self.nel_impact_hist_scrub)
                histf.compute_hist(x_impact, nel_impact * E_impact_eV,
                                   bias_x_hist, Dx_hist, self.energ_eV_impact_hist)

                # angle histogram
                if self.flag_cos_angle_hist:
                    histf.compute_hist(
                        costheta_impact, nel_impact, 0., self.cos_angle_width, self.cos_angle_hist)

                if flag_seg:
                    segi.update_seg_impact(
                        i_found, nel_impact, self.nel_hist_impact_seg)
                    segi.update_seg_impact(
                        i_found, nel_impact * E_impact_eV, self.energ_eV_impact_seg)

                    if self.flag_En_hist_seg:
                        for iseg in range(self.chamb.N_vert):
                            mask_this_seg = i_found == iseg
                            if np.sum(mask_this_seg) > 0:
                                En_imp_hist_this_seg = E_impact_eV[mask_this_seg]
                                En_imp_hist_this_seg[En_imp_hist_this_seg >
                                                     En_hist_max] = En_hist_max
                                histf.compute_hist(En_imp_hist_this_seg, nel_impact[mask_this_seg], 0., DEn_hist,
                                                   self.seg_En_hist_lines[iseg])

                En_imp_hist = E_impact_eV.copy()
                En_imp_hist[En_imp_hist > En_hist_max] = En_hist_max
                histf.compute_hist(En_imp_hist, nel_impact,
                                   0., DEn_hist, self.En_hist_line)

                self.Nel_impact_last_step = np.sum(nel_impact)
                self.En_imp_last_step_eV = np.sum(E_impact_eV * nel_impact)

                # Call secondary emission model
                (nel_emit_tot_events, event_type, event_info,
                    nel_replace, x_replace, y_replace, z_replace,
                    vx_replace, vy_replace, vz_replace, i_seg_replace,
                    nel_new_MPs, x_new_MPs, y_new_MPs, z_new_MPs,
                    vx_new_MPs, vy_new_MPs, vz_new_MPs, i_seg_new_MPs,
                 ) = sey_mod.impacts_on_surface(
                    MP_e.mass, nel_impact, x_impact, y_impact, z_impact,
                    vx_impact, vy_impact, vz_impact, Norm_x, Norm_y, i_found,
                    v_impact_n, E_impact_eV, costheta_impact, nel_mp_th, flag_seg
                )

                self.Nel_emit_last_step = np.sum(nel_emit_tot_events)

                # Replace old MPs
                x_mp[idx] = x_replace
                y_mp[idx] = y_replace
                z_mp[idx] = z_replace
                vx_mp[idx] = vx_replace
                vy_mp[idx] = vy_replace
                vz_mp[idx] = vz_replace
                nel_mp[idx] = nel_replace

                # subtract replaced macroparticles
                v_replace_mod = np.sqrt(
                    vx_replace**2 + vy_replace**2 + vz_replace**2)
                E_replace_eV = 0.5 * MP_e.mass / qe * v_replace_mod * v_replace_mod

                self.En_emit_last_step_eV = np.sum(E_replace_eV * nel_replace)

                histf.compute_hist(x_replace, -nel_replace * E_replace_eV,
                                   bias_x_hist, Dx_hist, self.energ_eV_impact_hist)
                if flag_seg:
                    segi.update_seg_impact(
                        i_seg_replace, -nel_replace * E_replace_eV, self.energ_eV_impact_seg)
                    segi.update_seg_impact(
                        i_seg_replace, nel_replace, self.nel_hist_emit_seg)

                # New macroparticles
                N_new_MPs = len(nel_new_MPs)
                if N_new_MPs > 0:
                    MP_e.add_new_MPs(N_new_MPs, nel_new_MPs, x_new_MPs, y_new_MPs, z_new_MPs,
                                     vx_new_MPs, vy_new_MPs, vz_new_MPs, tt_curr)

                    # subtract new macroparticles
                    v_new_MPs_mod = np.sqrt(
                        vx_new_MPs**2 + vy_new_MPs**2 + vz_new_MPs**2)
                    E_new_MPs_eV = 0.5 * MP_e.mass / qe * v_new_MPs_mod * v_new_MPs_mod

                    histf.compute_hist(x_new_MPs, -nel_new_MPs * E_new_MPs_eV,
                                       bias_x_hist, Dx_hist, self.energ_eV_impact_hist)

                    if flag_seg:
                        segi.update_seg_impact(
                            i_seg_new_MPs, -nel_new_MPs * E_new_MPs_eV, self.energ_eV_impact_seg)
                        segi.update_seg_impact(
                            i_seg_new_MPs, nel_new_MPs, self.nel_hist_emit_seg)

                    self.En_emit_last_step_eV += np.sum(
                        E_new_MPs_eV * nel_new_MPs)

        return MP_e

    def extract_sey_curves(self, n_rep, E_impact_eV_test, cos_theta_test, charge, mass):

        deltas = {}
        for etype in list(self.sey_mod.event_types.keys()):
            etype_name = self.sey_mod.event_types[etype]
            deltas[etype_name] = np.zeros(
                (len(cos_theta_test), len(E_impact_eV_test)))
        print('Extracting SEY curves...')
        for i_ct, ct in enumerate(cos_theta_test):
            print(('%d/%d' % (i_ct + 1, len(cos_theta_test))))
            for i_ene, Ene in enumerate(E_impact_eV_test):

                # nel_emit, flag_elast, flag_truesec = sey_mod.SEY_process(nel_impact=np.ones(n_rep),
                #                 E_impact_eV=Ene*np.ones(n_rep), costheta_impact=np.ones(n_rep)*ct, i_impact=np.array(n_rep*[0]))
                nel_impact = np.ones(n_rep)
                # Assuming normal is along x
                v_mod = np.sqrt(2 * Ene * qe / mass) * np.ones_like(nel_impact)
                vx = v_mod * ct
                vy = v_mod * np.sqrt(1 - ct * ct)

                nel_emit_tot_events, event_type, event_info,\
                    nel_replace, x_replace, y_replace, z_replace, vx_replace, vy_replace, vz_replace, i_seg_replace,\
                    nel_new_MPs, x_new_MPs, y_new_MPs, z_new_MPs, vx_new_MPs, vy_new_MPs, vz_new_MPs, i_seg_new_MPs =\
                    self.sey_mod.impacts_on_surface(
                        mass=mass, nel_impact=nel_impact, x_impact=nel_impact * 0, y_impact=nel_impact * 0, z_impact=nel_impact * 0,
                        vx_impact=vx * np.ones_like(nel_impact),
                        vy_impact=vy * np.ones_like(nel_impact),
                        vz_impact=nel_impact * 0,
                        Norm_x=np.ones_like(nel_impact), Norm_y=np.zeros_like(nel_impact),
                        i_found=np.int_(np.ones_like(nel_impact)),
                        v_impact_n=vx * np.ones_like(nel_impact),
                        E_impact_eV=Ene * np.ones_like(nel_impact),
                        costheta_impact=ct * np.ones_like(nel_impact),
                        nel_mp_th=1,
                        flag_seg=True)

                for etype in list(self.sey_mod.event_types.keys()):
                    etype_name = self.sey_mod.event_types[etype]
                    thisdelta = deltas[etype_name]
                    thisdelta[i_ct, i_ene] = np.sum(
                        nel_emit_tot_events[event_type == etype]) / np.sum(nel_impact)
                    deltas[etype_name] = thisdelta

        print('Done extracting SEY curves.')

        return deltas

    def extract_energy_distributions(self, n_rep, E_impact_eV_test, cos_theta_test, mass, Nbin_extract_ene, factor_ene_dist_max):
        """Extract energy distributions for secondary electrons."""
        emit_ene_g_hist = np.linspace(
            0., E_impact_eV_test * factor_ene_dist_max, Nbin_extract_ene)
        Dextract_ene = emit_ene_g_hist[1] - emit_ene_g_hist[0]
        extract_ene_hist = {}

        for etype in list(self.sey_mod.event_types.keys()):
            etype_name = self.sey_mod.event_types[etype]
            extract_ene_hist[etype_name] = np.zeros(
                shape=(len(emit_ene_g_hist), len(cos_theta_test)), dtype=float)

        print('Extracting energy distributions...')
        for i_ct, ct in enumerate(cos_theta_test):
            print(('%d/%d' % (i_ct + 1, len(cos_theta_test))))
            Ene = E_impact_eV_test
            nel_impact = np.ones(n_rep)
            # Assuming normal is along x
            v_mod = np.sqrt(2 * Ene * qe / mass) * np.ones_like(nel_impact)
            vx = v_mod * ct
            vy = v_mod * np.sqrt(1 - ct * ct)

            nel_emit_tot_events, event_type, event_info,\
                nel_replace, x_replace, y_replace, z_replace, vx_replace, vy_replace, vz_replace, i_seg_replace,\
                nel_new_MPs, x_new_MPs, y_new_MPs, z_new_MPs, vx_new_MPs, vy_new_MPs, vz_new_MPs, i_seg_new_MPs =\
                self.sey_mod.impacts_on_surface(
                    mass=mass, nel_impact=nel_impact, x_impact=nel_impact * 0, y_impact=nel_impact * 0, z_impact=nel_impact * 0,
                    vx_impact=vx * np.ones_like(nel_impact),
                    vy_impact=vy * np.ones_like(nel_impact),
                    vz_impact=nel_impact * 0,
                    Norm_x=np.ones_like(nel_impact), Norm_y=np.zeros_like(nel_impact),
                    i_found=np.int_(np.ones_like(nel_impact)),
                    v_impact_n=vx * np.ones_like(nel_impact),
                    E_impact_eV=Ene * np.ones_like(nel_impact),
                    costheta_impact=ct * np.ones_like(nel_impact),
                    nel_mp_th=1,
                    flag_seg=True)

            v_replace_mod = np.sqrt(
                vx_replace**2 + vy_replace**2 + vz_replace**2)
            E_replace_eV = 0.5 * mass / qe * v_replace_mod * v_replace_mod

            v_new_MPs_mod = np.sqrt(
                vx_new_MPs**2 + vy_new_MPs**2 + vz_new_MPs**2)
            E_new_MPs_eV = 0.5 * mass / qe * v_new_MPs_mod * v_new_MPs_mod

            E_all_MPs_eV = np.concatenate([E_replace_eV, E_new_MPs_eV])

            extended_event_type = event_info['extended_event_type']
            for etype in list(self.sey_mod.event_types.keys()):
                etype_name = self.sey_mod.event_types[etype]
                extract_type = extract_ene_hist[etype_name]
                # if there are no events of type etype
                if E_all_MPs_eV[extended_event_type == etype].shape == (0,):
                    pass
                else:
                    temp = extract_type[:, i_ct].copy()
                    histf.compute_hist(E_all_MPs_eV[extended_event_type == etype], np.ones(
                        len(E_all_MPs_eV[extended_event_type == etype])), 0., Dextract_ene, temp)
                    extract_type[:, i_ct] = temp
                extract_ene_hist[etype_name] = extract_type

        extract_ene_hist['emit_ene_g_hist'] = emit_ene_g_hist

        print('Done extracting energy distributions.')

        return extract_ene_hist

#----------------------------------------------------------------------
#----------------------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 6.7.0
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

import numpy as np
from scipy.constants import c, e
import time
from PyEC4PyHT import Ecloud
from gas_ionization_class import residual_gas_ionization


class MP_light(object):
    pass



class Ecloud_fastion(Ecloud):

    def __init__(self, L_ecloud, slicer, Dt_ref, pyecl_input_folder = './',
                flag_clean_slices = False, slice_by_slice_mode = False, space_charge_obj = None,
                beam_monitor = None, include_cloud_sc = False, ionize_only_first_bunch = False, **kwargs):


        super(Ecloud_fastion, self).__init__(L_ecloud, slicer, Dt_ref, pyecl_input_folder = pyecl_input_folder,
                                                flag_clean_slices = flag_clean_slices, slice_by_slice_mode = slice_by_slice_mode,
                                                space_charge_obj = space_charge_obj, **kwargs)

        self.beam_monitor = beam_monitor
        self.include_cloud_sc = include_cloud_sc
        self.ionize_only_first_bunch = ionize_only_first_bunch

        self.MP_e_field_state = self.spacech_ele.PyPICobj.get_state_object()
        self.MP_p_field_state = self.spacech_ele.PyPICobj.get_state_object()

        self.gas_ionization = self.resgasion



    #@profile
    def track(self, beam):

        start_time = time.mktime(time.localtime())

        if self.track_only_first_time:
            if self.N_tracks>0:
                print 'Warning: Track skipped because track_only_first_time is True.'
                return

        self._reinitialize()

        if hasattr(beam.particlenumber_per_mp, '__iter__'):
            raise ValueError('ecloud module assumes same size for all beam MPs')

        if self.flag_clean_slices:
            beam.clean_slices()

        slices = beam.get_slices(self.slicer)
        self.slicer.add_statistics(sliceset=slices, beam=beam, statistics=True)

        # Only track over slices with particles
        filled_slices = np.where(slices.n_macroparticles_per_slice > 0)[0]

        for i in filled_slices[::-1]:

            # select particles in the bunch
            ix = slices.particle_indices_of_slice(i)

            # slice size and time step
            dz = (slices.z_bins[i + 1] - slices.z_bins[i]) # in this case it is the bucket length

            self._track_single_slice(beam, ix, dz)

        if self.beam_monitor is not None:
            self.beam_monitor.dump(beam)

        self._finalize()

        self.N_tracks+=1

        stop_time = time.mktime(time.localtime())
        print 'Done track in ', (stop_time-start_time), 's'


    #@profile
    def _track_in_single_slice_mode(self, beam):

        if hasattr(beam.particlenumber_per_mp, '__iter__'):
            raise ValueError('ecloud module assumes same size for all beam MPs')

        if self.flag_clean_slices:
            raise ValueError(
                    'track cannot clean the slices in slice-by-slice mode! ')

        if beam.slice_info is not 'unsliced': # and beam.macroparticlenumber > 0:
            dz = beam.slice_info['z_bin_right']-beam.slice_info['z_bin_left']
            self._track_single_slice(beam, ix=np.arange(beam.macroparticlenumber), dz=dz)



    def generate_twin_ecloud_with_shared_space_charge(self):
        if hasattr(self, 'efieldmap'):
            raise ValueError('Ecloud has been replaced with field map. I cannot generate a twin ecloud!')

        return Ecloud_fastion(self.L_ecloud, self.slicer, self.Dt_ref, self.MP_e.mass, self.MP_e.charge,
                self.pyecl_input_folder, flag_clean_slices = self.flag_clean_slices,
                slice_by_slice_mode = self.slice_by_slice_mode, space_charge_obj = self.spacech_ele,
                beam_monitor = self.beam_monitor, include_cloud_sc = self.include_cloud_sc,
                ionize_only_first_bunch = self.ionize_only_first_bunch, **self.kwargs)


    #@profile
    def _track_single_slice(self, beam, ix, dz):

        #pass
        if len(ix) > 0:

            MP_e = self.MP_e
            dynamics = self.dynamics
            impact_man = self.impact_man
            spacech_ele = self.spacech_ele
            MP_e_state = self.MP_e_field_state
            MP_p_state = self.MP_p_field_state

            dt = dz / (beam.beta * c)

            # define substep
            if dt > self.Dt_ref:
                N_sub_steps = int(np.round(dt / self.Dt_ref))
            else:
                N_sub_steps = 1

            Dt_substep = dt / N_sub_steps
            # print Dt_substep, N_sub_steps, dt

            # beam particles
            MP_p = MP_light()
            MP_p.x_mp = beam.x[ix]
            MP_p.y_mp = beam.y[ix]
            MP_p.nel_mp = beam.x[ix] * 0. + beam.particlenumber_per_mp
            MP_p.N_mp = len(beam.x[ix])
            MP_p.charge = beam.charge

            mean_x = np.mean(beam.x[ix])
            mean_y = np.mean(beam.y[ix])
            sigma_x = np.std(beam.x[ix])
            sigma_y = np.std(beam.y[ix])

            if self.gas_ion_flag == 1:
                Np_bunch = MP_p.N_mp * beam.particlenumber_per_mp
                dz_bunch = dz
                lambda_bunch = Np_bunch
                dt_bunch = 1 / c
                MP_e = self.gas_ionization.generate(MP_e=MP_e, lambda_t=lambda_bunch, Dt=dt_bunch, sigmax=sigma_x,
                                                    sigmay=sigma_y, x_beam_pos=mean_x, y_beam_pos=mean_y)
                if self.ionize_only_first_bunch:
                    self.gas_ion_flag = 0

            # scatter fields
            MP_e_state.scatter(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp],MP_e.nel_mp[0:MP_e.N_mp], charge = MP_e.charge)
            MP_p_state.scatter(MP_p.x_mp[0:MP_p.N_mp],MP_p.y_mp[0:MP_p.N_mp],MP_p.nel_mp[0:MP_p.N_mp], charge = MP_p.charge)

            # solve fields
            spacech_ele.PyPICobj.solve_states([MP_e_state, MP_p_state])

            # gather fields
            Ex_sc_p, Ey_sc_p = MP_e_state.gather(MP_p.x_mp[0:MP_p.N_mp],MP_p.y_mp[0:MP_p.N_mp])
            Ex_n_beam, Ey_n_beam = MP_p_state.gather(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp])

            # kick cloud particles
            MP_e.vx_mp[:MP_e.N_mp] += Ex_n_beam * MP_e.charge / MP_e.mass / c
            MP_e.vy_mp[:MP_e.N_mp] += Ey_n_beam * MP_e.charge / MP_e.mass / c

            # kick beam particles
            fact_kick = beam.charge / (beam.mass * beam.beta * beam.beta * beam.gamma * c * c) * self.L_ecloud
            beam.xp[ix] += fact_kick * Ex_sc_p
            beam.yp[ix] += fact_kick * Ey_sc_p


            # Total electric field on electrons
            if self.include_cloud_sc:
                Ex_sc_n, Ey_sc_n = MP_e_state.gather(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp])
                Ex_n = Ex_sc_n
                Ey_n = Ey_sc_n
            else:
                Ex_n = MP_e.vx_mp[:MP_e.N_mp] * 0.
                Ey_n = Ex_n

            # save position before motion step
            old_pos = MP_e.get_positions()

            # motion electrons
            MP_e = dynamics.stepcustomDt(MP_e, Ex_n,Ey_n, Dt_substep=Dt_substep, N_sub_steps=N_sub_steps)

            # impacts: backtracking and secondary emission
            MP_e = impact_man.backtrack_and_second_emiss(old_pos, MP_e)


            if self.save_ele_distributions_last_track:
                self.rho_ele_last_track.append(spacech_ele.rho.copy())
                #print 'Here'

            if self.save_ele_potential_and_field:
                self.phi_ele_last_track.append(spacech_ele.phi.copy())
                self.Ex_ele_last_track.append(spacech_ele.efx.copy())
                self.Ey_ele_last_track.append(spacech_ele.efy.copy())

            if self.save_ele_MP_position:
                self.x_MP_last_track.append(MP_e.x_mp.copy())
                self.y_MP_last_track.append(MP_e.y_mp.copy())

            if self.save_ele_MP_velocity:
                self.vx_MP_last_track.append(MP_e.vx_mp.copy())
                self.vy_MP_last_track.append(MP_e.vy_mp.copy())

            if self.save_ele_MP_size:
                self.nel_MP_last_track.append(MP_e.nel_mp.copy())

            if self.save_ele_MP_position or self.save_ele_MP_velocity or self.save_ele_MP_size:
                self.N_MP_last_track.append(MP_e.N_mp)

        else:
            pass

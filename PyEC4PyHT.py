#-Begin-preamble-------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 7.4.0
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

import os
import subprocess
import time

import numpy as np
from scipy.constants import c, e, m_e

import myloadmat_to_obj as mlm
import init
import buildup_simulation as bsim

class Empty(object):
    pass

class DummyBeamTim(object):
    def __init__(self, PyPIC_state):
        self.PyPIC_state = PyPIC_state

        self.b_spac = 0.
        self.pass_numb = 0
        self.N_pass_tot = 1

    def get_beam_eletric_field(self, MP_e):
        if (MP_e.N_mp>0):
            if self.PyPIC_state is None:
                Ex_n_beam = 0. * MP_e.x_mp[0:MP_e.N_mp]
                Ey_n_beam = 0. * MP_e.y_mp[0:MP_e.N_mp]
            else:
                # compute beam electric field
                Ex_n_beam, Ey_n_beam = self.PyPIC_state.gather(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp])
        else:
            Ex_n_beam=0.
            Ey_n_beam=0.
        return Ex_n_beam, Ey_n_beam


extra_allowed_kwargs = {'x_beam_offset', 'y_beam_offset', 'probes_position'}

class Ecloud(object):
    
    
    def generate_twin_ecloud_with_shared_space_charge(self):
        if hasattr(self, 'efieldmap'):
            raise ValueError('Ecloud has been replaced with field map. I cannot generate a twin ecloud!')
        return Ecloud(self.L_ecloud, self.slicer, self.Dt_ref, self.pyecl_input_folder, self.flag_clean_slices,
                self.slice_by_slice_mode, self.spacech_ele, self.kick_mode_for_beam_field,
                self.beam_monitor, self.verbose, self.save_pyecl_outp_as, **self.kwargs)

    
    def __init__(self, L_ecloud, slicer, Dt_ref, pyecl_input_folder='./', flag_clean_slices=False,
                 slice_by_slice_mode=False, space_charge_obj=None, kick_mode_for_beam_field=False,
                 beam_monitor=None, verbose=False, save_pyecl_outp_as=None, 
                 **kwargs):

        print 'PyECLOUD Version 7.4.0'

        # These git commands return the hash and the branch of the specified git directory.
        path_to_git = os.path.dirname(os.path.abspath(__file__)) +'/.git'
        cmd_hash = 'git --git-dir %s rev-parse HEAD' % path_to_git
        cmd_branch = 'git --git-dir %s rev-parse --abbrev-ref HEAD' % path_to_git
        try:
            git_hash = 'git hash: %s' % (subprocess.check_output(cmd_hash.split()).split()[0])
        except Exception as e:
            git_hash = 'Retrieving git hash failed'
            print(e)
        print(git_hash)

        try:
            git_branch = 'git branch: %s' % (subprocess.check_output(cmd_branch.split()).split()[0])
        except Exception as e:
            git_branch = 'Retrieving git branch failed'
            print(e)
        print(git_branch)

        print 'PyHEADTAIL module'
        print 'Initializing ecloud from folder: '+pyecl_input_folder
        self.slicer = slicer
        self.Dt_ref = Dt_ref
        self.L_ecloud = L_ecloud

        self.pyecl_input_folder = pyecl_input_folder
        self.kwargs = kwargs

        self.cloudsim = bsim.BuildupSimulation(
                    pyecl_input_folder=pyecl_input_folder, skip_beam=True, 
                    spacech_ele=space_charge_obj,
                    ignore_kwargs=extra_allowed_kwargs, 
                    skip_pyeclsaver=(save_pyecl_outp_as is None), 
                    filen_main_outp = save_pyecl_outp_as,
                    **self.kwargs)


        if self.cloudsim.config_dict['track_method'] == 'Boris':
            pass
        elif self.cloudsim.config_dict['track_method'] == 'BorisMultipole':
            pass
        else:
            raise ValueError("""track_method should be 'Boris' or 'BorisMultipole' - others are not implemented in the PyEC4PyHT module""")

        self.x_beam_offset = 0.
        self.y_beam_offset = 0.
        if 'x_beam_offset' in kwargs:
            self.x_beam_offset = kwargs['x_beam_offset']
        if 'y_beam_offset' in kwargs:
            self.y_beam_offset = kwargs['y_beam_offset']


        # initialize proton density probes
        self.save_ele_field_probes = False
        self.x_probes = -1
        self.y_probes = -1
        self.Ex_ele_last_track_at_probes = -1
        self.Ey_ele_last_track_at_probes = -1
        if 'probes_position' in kwargs.keys():
            self.save_ele_field_probes = True
            self.probes_position = kwargs['probes_position']
            self.N_probes = len(self.probes_position)
            self.x_probes = []
            self.y_probes = []
            for ii_probe in xrange(self.N_probes):
                self.x_probes.append(self.probes_position[ii_probe]['x'])
                self.y_probes.append(self.probes_position[ii_probe]['y'])

            self.x_probes = np.array(self.x_probes)
            self.y_probes = np.array(self.y_probes)

        self.N_tracks = 0

        #self.cloudsim.spacech_ele.flag_decimate = False



        if self.cloudsim.flag_multiple_clouds:
            raise ValueError('Multiple clouds not yet tested in PyEC4PyHT!')


        self.save_ele_distributions_last_track = False
        self.save_ele_potential_and_field = False
        self.save_ele_potential = False
        self.save_ele_field = False
        self.save_ele_MP_position = False
        self.save_ele_MP_velocity = False
        self.save_ele_MP_size = False

        self.track_only_first_time = False

        self.initial_MP_e_clouds = [cl.MP_e.extract_dict() for cl in self.cloudsim.cloud_list]

        self.flag_clean_slices = flag_clean_slices

        self.beam_PyPIC_state = self.cloudsim.spacech_ele.PyPICobj.get_state_object()

        self.slice_by_slice_mode = slice_by_slice_mode
        if self.slice_by_slice_mode:
            self.track = self._track_in_single_slice_mode
            self.finalize_and_reinitialize = self._finalize_and_reinitialize

        self.spacech_ele = self.cloudsim.spacech_ele # For backwards compatibility

        self.kick_mode_for_beam_field = kick_mode_for_beam_field
        self.beam_monitor = beam_monitor
        
        self.verbose = verbose
        self.save_pyecl_outp_as = save_pyecl_outp_as

        self.i_reinit = 0
        self.t_sim = 0.
        self.i_curr_bunch = -1



    #    @profile
    def track(self, beam):

        if self.track_only_first_time:
            if self.N_tracks>0:
                print 'Warning: Track skipped because track_only_first_time is True.'
                return
        
        if self.verbose:
            start_time = time.mktime(time.localtime())
        
        self._reinitialize()

        if hasattr(beam.particlenumber_per_mp, '__iter__'):
            raise ValueError('ecloud module assumes same size for all beam MPs')

        if self.flag_clean_slices:
            beam.clean_slices()

        slices = beam.get_slices(self.slicer)

        for i in xrange(slices.n_slices-1, -1, -1):

            # select particles in the slice
            ix = slices.particle_indices_of_slice(i)

            # slice size and time step
            dz = (slices.z_bins[i + 1] - slices.z_bins[i])

            self._track_single_slice(beam, ix, dz, force_pyecl_newpass=(i==0))
            
        
        # Used by Lotta to debug fastion mode
        if self.beam_monitor is not None:
            self.beam_monitor.dump(beam)

        self._finalize()
        
        if self.verbose:
            stop_time = time.mktime(time.localtime())
            print 'Done track %d in %.1f s'%(self.N_tracks, stop_time-start_time)

        self.N_tracks+=1

    def replace_with_recorded_field_map(self, delete_ecloud_data=True):

        if self.track_only_first_time:
            print 'Warning: replace_with_recorded_field_map resets track_only_first_time = False'
            self.track_only_first_time=False

        if not hasattr(self, 'efieldmap'):
            from Transverse_Efield_map_for_frozen_cloud import Transverse_Efield_map
            self.efieldmap = Transverse_Efield_map(xg = self.spacech_ele.xg, yg = self.spacech_ele.yg,
                Ex=self.Ex_ele_last_track, Ey=self.Ey_ele_last_track, L_interaction=self.L_ecloud,
                slicer = self.slicer,
                flag_clean_slices = True,
                x_beam_offset = self.x_beam_offset, y_beam_offset = self.y_beam_offset,
                slice_by_slice_mode = self.slice_by_slice_mode)

            self._ecloud_track = self.track

            self.track = self.efieldmap.track
            self.finalize_and_reinitialize = self.efieldmap.finalize_and_reinitialize

            if delete_ecloud_data:
                self.spacech_ele=None
                self.cloudsim=None

        else:
            print 'Warning: efieldmap already exists. I do nothing.'

    def track_once_and_replace_with_recorded_field_map(self, bunch, delete_ecloud_data=True):
        self.save_ele_field = True
        self.track_only_first_time = True
        if self.slice_by_slice_mode:
            if not hasattr(bunch, '__iter__'):
                raise ValueError('A list of slices should be provided!')
            self._reinitialize()
            for slc in bunch:
                self.track(slc)
            self._finalize()
        else:
            self.track(bunch)
        self.save_ele_field = False
        self.track_only_first_time = False
        self.replace_with_recorded_field_map(delete_ecloud_data=delete_ecloud_data)


    def _track_single_slice(self, slic, ix, dz, force_pyecl_newpass=False):

        spacech_ele = self.cloudsim.spacech_ele

        # Check if the slice interacts with the beam
        if hasattr(slic, 'slice_info'):
            if 'interact_with_EC' in slic.slice_info.keys():
                interact_with_EC = slic.slice_info['interact_with_EC']
            else:
                interact_with_EC = True
        else:
            interact_with_EC = True

        # Compute slice length
        dt_slice = dz / (slic.beta * c)

        # Check if sub-slicing is needed
        if self.cloudsim.config_dict['Dt'] is not None:
            if dt_slice>self.cloudsim.config_dict['Dt']:
                if interact_with_EC: 
                    raise ValueError('Slices that interact with the cloud cannot be longer than the buildup timestep!')

                N_cloud_steps = np.int_(np.ceil(dt_slice/self.cloudsim.config_dict['Dt']))
                dt_cloud_step = dt_slice/N_cloud_steps
                dt_array = np.array(N_cloud_steps*[dt_cloud_step])
            else:
                dt_array = np.array([dt_slice])
        else:
            dt_array = np.array([dt_slice])

        # Acquire bunch passage information
        if hasattr(slic, 'slice_info'):
            if 'info_parent_bunch' in slic.slice_info.keys():
                # check if new passage
                if slic.slice_info['info_parent_bunch']['i_bunch']>self.i_curr_bunch:
                    self.i_curr_bunch = slic.slice_info['info_parent_bunch']['i_bunch']
                    new_pass = True
                else:
                    new_pass = force_pyecl_newpass
                    
                # check if first slice of first bunch
                if slic.slice_info['info_parent_bunch']['i_bunch']==0 and slic.slice_info['i_slice']==0:
                    self.finalize_and_reinitialize()                       
            else:
                new_pass = force_pyecl_newpass
                self.i_curr_bunch = 0
        else:
            new_pass = force_pyecl_newpass
            self.i_curr_bunch = 0


        # Cloud simulation
        for i_clou_step, dt in enumerate(dt_array): 

            # define substep
            if dt>self.Dt_ref:
                N_sub_steps = int(np.round(dt/self.Dt_ref))
            else:
                N_sub_steps=1

            Dt_substep = dt/N_sub_steps
            #print Dt_substep, N_sub_steps, dt

            if len(ix)==0 or not interact_with_EC: #no particles in the beam
                
                #build dummy beamtim object
                dummybeamtim = DummyBeamTim(None)

                dummybeamtim.lam_t_curr = 0.
                dummybeamtim.sigmax = 0. 
                dummybeamtim.sigmay = 0.
                dummybeamtim.x_beam_pos = 0.
                dummybeamtim.y_beam_pos = 0.
            else:
                
                # beam field
                self.beam_PyPIC_state.scatter(
                            x_mp = slic.x[ix]+self.x_beam_offset, 
                            y_mp = slic.y[ix]+self.y_beam_offset, 
                            nel_mp = slic.x[ix]*0.+slic.particlenumber_per_mp/dz,
                            charge = slic.charge)
                self.cloudsim.spacech_ele.PyPICobj.solve_states([self.beam_PyPIC_state])

                #build dummy beamtim object
                dummybeamtim = DummyBeamTim(self.beam_PyPIC_state)

                dummybeamtim.lam_t_curr = np.mean(slic.particlenumber_per_mp/dz)*len(ix)
                dummybeamtim.sigmax = np.std(slic.x[ix]) 
                dummybeamtim.sigmay = np.std(slic.y[ix])
                dummybeamtim.x_beam_pos = np.mean(slic.x[ix])+self.x_beam_offset
                dummybeamtim.y_beam_pos = np.mean(slic.y[ix])+self.y_beam_offset 
            
            dummybeamtim.tt_curr = self.t_sim # In order to have the PIC activated
            dummybeamtim.Dt = dt
            dummybeamtim.pass_numb = self.i_curr_bunch
            dummybeamtim.flag_new_bunch_pass = new_pass

            # Force space charge recomputation 
            force_recompute_space_charge = interact_with_EC or (i_clou_step==0) # we always force at the first step, as we don't know the state of the PIC
            
            # Disable cleanings and regenerations
            skip_MP_cleaning = interact_with_EC
            skip_MP_regen = interact_with_EC

            # print(dummybeamtim.tt_curr, dummybeamtim.flag_new_bunch_pass, force_recompute_space_charge)

            # Perform cloud simulation step
            self.cloudsim.sim_time_step(beamtim_obj=dummybeamtim, 
                    Dt_substep_custom=Dt_substep, N_sub_steps_custom=N_sub_steps, 
                    kick_mode_for_beam_field=self.kick_mode_for_beam_field,
                    force_recompute_space_charge=force_recompute_space_charge,
                    skip_MP_cleaning=skip_MP_cleaning, skip_MP_regen=skip_MP_regen)


            # Build MP_system-like object with beam coordinates
            MP_p = Empty()
            MP_p.x_mp = slic.x[ix]+self.x_beam_offset
            MP_p.y_mp = slic.y[ix]+self.y_beam_offset
            MP_p.N_mp = len(slic.x[ix])

            ## compute cloud field on beam particles
            Ex_sc_p, Ey_sc_p = spacech_ele.get_sc_eletric_field(MP_p)

            ## kick beam particles
            fact_kick = slic.charge/(slic.mass*slic.beta*slic.beta*slic.gamma*c*c)*self.L_ecloud
            slic.xp[ix]+=fact_kick*Ex_sc_p
            slic.yp[ix]+=fact_kick*Ey_sc_p

            ## Diagnostics
            if self.save_ele_distributions_last_track:
                self.rho_ele_last_track.append(spacech_ele.rho.copy())
                #print 'Here'

            if self.save_ele_potential:
                self.phi_ele_last_track.append(spacech_ele.phi.copy())

            if self.save_ele_field:
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

            if self.save_ele_field_probes:
                MP_probes = Empty()
                MP_probes.x_mp = self.x_probes
                MP_probes.y_mp = self.y_probes
                MP_probes.nel_mp = self.x_probes*0.+1. #fictitious charge of 1 C
                MP_probes.N_mp = len(self.x_probes)
                Ex_sc_probe, Ey_sc_probe = spacech_ele.get_sc_eletric_field(MP_probes)

                self.Ex_ele_last_track_at_probes.append(Ex_sc_probe.copy())
                self.Ey_ele_last_track_at_probes.append(Ey_sc_probe.copy())

            self.t_sim += dt
            new_pass = False # it can be true only for the first sub-slice of a slice

    def _reinitialize(self):

        cc = mlm.obj_from_dict(self.cloudsim.config_dict)
        
        
        for thiscloud, initdict in zip(self.cloudsim.cloud_list, self.initial_MP_e_clouds):
            thiscloud.MP_e.init_from_dict(initdict)
            thiscloud.MP_e.set_nel_mp_ref(thiscloud.MP_e.nel_mp_ref_0)
            
            
            thisconf = thiscloud.config_dict
            
            if thiscloud.pyeclsaver is not None:
                thiscloud.pyeclsaver.start_observing(cc.Dt, thiscloud.MP_e, None, thiscloud.impact_man,
                                       thisconf.r_center, thisconf.Dt_En_hist, thisconf.logfile_path, thisconf.progress_path, 
                                       flag_detailed_MP_info=thisconf.flag_detailed_MP_info, flag_movie=thisconf.flag_movie, 
                                       flag_sc_movie=thisconf.flag_sc_movie, save_mp_state_time_file=thisconf.save_mp_state_time_file,
                                       flag_presence_sec_beams=self.cloudsim.flag_presence_sec_beams, sec_beams_list=self.cloudsim.sec_beams_list, dec_fac_secbeam_prof=thisconf.dec_fac_secbeam_prof,
                                       el_density_probes=thisconf.el_density_probes, save_simulation_state_time_file=thisconf.save_simulation_state_time_file,
                                       x_min_hist_det=thisconf.x_min_hist_det, x_max_hist_det=thisconf.x_max_hist_det, 
                                       y_min_hist_det=thisconf.y_min_hist_det, y_max_hist_det=thisconf.y_max_hist_det,
                                       Dx_hist_det=thisconf.Dx_hist_det, dec_fact_out=cc.dec_fact_out, stopfile=cc.stopfile, filen_main_outp=thisconf.filen_main_outp,
                                       flag_cos_angle_hist=thisconf.flag_cos_angle_hist, cos_angle_width=thisconf.cos_angle_width,
                                       flag_multiple_clouds=(len(self.cloudsim.cloud_list)>1), cloud_name=thisconf.cloud_name, flag_last_cloud=(thiscloud is self.cloudsim.cloud_list[-1]))
            
                thiscloud.pyeclsaver.filen_main_outp = thiscloud.pyeclsaver.filen_main_outp.split('.mat')[0].split('__iter')[0]+'__iter%d.mat'%self.i_reinit
            

        if self.save_ele_distributions_last_track:
            self.rho_ele_last_track = []

        if self.save_ele_potential_and_field:
            self.save_ele_potential = True
            self.save_ele_field = True

        if self.save_ele_potential:
            self.phi_ele_last_track = []

        if self.save_ele_field:
            self.Ex_ele_last_track = []
            self.Ey_ele_last_track = []

        if self.save_ele_MP_position:
            self.x_MP_last_track = []
            self.y_MP_last_track = []

        if self.save_ele_MP_velocity:
            self.vx_MP_last_track = []
            self.vy_MP_last_track = []

        if self.save_ele_MP_size:
            self.nel_MP_last_track = []

        if self.save_ele_MP_position or self.save_ele_MP_velocity or self.save_ele_MP_size:
            self.N_MP_last_track = []

        if self.save_ele_field_probes:
            self.Ex_ele_last_track_at_probes = []
            self.Ey_ele_last_track_at_probes = []

        #~ self.t_sim = 0.
        self.i_curr_bunch = -1
        self.i_reinit += 1


    def _finalize(self):

        if self.save_ele_distributions_last_track:
            self.rho_ele_last_track = np.array(self.rho_ele_last_track[::-1])

        if self.save_ele_potential:
            self.phi_ele_last_track = np.array(self.phi_ele_last_track[::-1])

        if self.save_ele_field:
            self.Ex_ele_last_track = np.array(self.Ex_ele_last_track[::-1])
            self.Ey_ele_last_track = np.array(self.Ey_ele_last_track[::-1])

        if self.save_ele_MP_position:
            self.x_MP_last_track = np.array(self.x_MP_last_track[::-1])
            self.y_MP_last_track = np.array(self.y_MP_last_track[::-1])

        if self.save_ele_MP_velocity:
            self.vx_MP_last_track = np.array(self.vx_MP_last_track[::-1])
            self.vy_MP_last_track = np.array(self.vy_MP_last_track[::-1])

        if self.save_ele_MP_size:
            self.nel_MP_last_track = np.array(self.nel_MP_last_track[::-1])

        if self.save_ele_MP_position or self.save_ele_MP_velocity or self.save_ele_MP_size:
                self.N_MP_last_track = np.array(self.N_MP_last_track[::-1])

        if self.save_ele_field_probes:
            self.Ex_ele_last_track_at_probes = np.array(self.Ex_ele_last_track_at_probes[::-1])
            self.Ey_ele_last_track_at_probes = np.array(self.Ey_ele_last_track_at_probes[::-1])

    def _finalize_and_reinitialize(self):
        print('Exec. finalize and reinitialize')
        self._finalize()
        self._reinitialize()

    def _track_in_single_slice_mode(self, beam):

        if hasattr(beam.particlenumber_per_mp, '__iter__'):
            raise ValueError('ecloud module assumes same size for all beam MPs')

        if self.flag_clean_slices:
            raise ValueError(
                    'track cannot clean the slices in slice-by-slice mode! ')

        if beam.slice_info != 'unsliced':
            dz = beam.slice_info['z_bin_right']-beam.slice_info['z_bin_left']
            self._track_single_slice(beam, ix=np.arange(beam.macroparticlenumber), dz=dz)
            
    def remove_savers(self):
        for thiscloud in self.cloudsim.cloud_list:
            thiscloud.pyeclsaver = None
            


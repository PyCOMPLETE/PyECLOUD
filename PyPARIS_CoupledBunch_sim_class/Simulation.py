# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

import types
import os
import importlib

import numpy as np
from scipy.constants import c

import PyPARIS.communication_helpers as ch
import PyPARIS.share_segments as shs
import PyPARIS.slicing_tool as sl
import PyECLOUD.myfilemanager as mfm

from PyHEADTAIL.particles.slicing import UniformBinSlicer

from . import Save_Load_Status as SLS

import Simulation_parameters as pp

flag_aperture = True # never tested otherwise


class Simulation(object):
    def __init__(self):
        self.N_turns = pp.N_turns
        self.N_buffer_float_size = 10000000
        self.N_buffer_int_size = 200
        self.N_parellel_rings = pp.N_parellel_rings
        
        self.n_slices_per_bunch = pp.n_slices_per_bunch
        self.z_cut_slicing = pp.z_cut_slicing
        self.N_pieces_per_transfer = pp.N_pieces_per_transfer 
        self.verbose = pp.verbose
        self.mpi_verbose = pp.mpi_verbose
        self.enable_barriers = pp.enable_barriers
        self.save_beam_at_turns = pp.save_beam_at_turns


    def pre_init_master(self):
        # Manage multi-run operation
        SimSt = SLS.SimulationStatus(N_turns_per_run=self.N_turns,
                resubmit_command = pp.resubmit_command, N_turns_target=pp.N_turns_target,
                submission_system = pp.submission_system if hasattr(pp, 'submission_system') else None)
        SimSt.before_simulation()
        return [SimSt.to_string()]
        

    def init_all(self, from_master):
        
        print('Exec init...')

        # Manage multi-run operation
        self.SimSt = SLS.SimulationStatus(N_turns_per_run=self.N_turns,
                resubmit_command = pp.resubmit_command, N_turns_target=pp.N_turns_target,
                submission_system = pp.submission_system if hasattr(pp, 'submission_system') else None)        
        self.SimSt.from_string(from_master[0])

        machine_name_strings = pp.machine_class_path.split('.')
        machine_module_name = '.'.join(machine_name_strings[:-1])
        machine_class_name = machine_name_strings[-1]
        machine_module = importlib.import_module(machine_module_name)
        machine_class = getattr(machine_module, machine_class_name)
        self.machine = machine_class(n_segments=pp.n_segments, machine_configuration=pp.machine_configuration,
                        Qp_x=pp.Qp_x, Qp_y=pp.Qp_y,
                        octupole_knob=pp.octupole_knob)

        self.n_non_parallelizable = 1 #RF

        inj_optics = self.machine.transverse_map.get_injection_optics()
        sigma_x_smooth = np.sqrt(inj_optics['beta_x']*pp.epsn_x/self.machine.betagamma)
        sigma_y_smooth = np.sqrt(inj_optics['beta_y']*pp.epsn_y/self.machine.betagamma)

        if flag_aperture:
            # setup transverse losses (to "protect" the ecloud)
            import PyHEADTAIL.aperture.aperture as aperture
            apt_xy = aperture.EllipticalApertureXY(x_aper=pp.target_size_internal_grid_sigma*sigma_x_smooth, 
                                                   y_aper=pp.target_size_internal_grid_sigma*sigma_x_smooth)
            self.machine.one_turn_map.append(apt_xy)
            self.n_non_parallelizable +=1 

        if pp.enable_transverse_damper:
            # setup transverse damper
            from PyHEADTAIL.feedback.transverse_damper import TransverseDamper
            damper = TransverseDamper(dampingrate_x=pp.dampingrate_x, dampingrate_y=pp.dampingrate_y)
            self.machine.one_turn_map.append(damper)
            self.n_non_parallelizable +=1
            
        if pp.enable_ecloud:
            print('Build ecloud...')
            import PyECLOUD.PyEC4PyHT as PyEC4PyHT
            ecloud = PyEC4PyHT.Ecloud(
                    L_ecloud=pp.L_ecloud_tot/pp.n_segments, slicer=None, slice_by_slice_mode=True,
                    Dt_ref=pp.Dt_ref, pyecl_input_folder=pp.pyecl_input_folder,
                    chamb_type=pp.chamb_type,
                    filename_chm=pp.filename_chm, 
                    target_grid={'x_min_target':-pp.target_size_internal_grid_sigma*sigma_x_smooth,
                                 'x_max_target':pp.target_size_internal_grid_sigma*sigma_x_smooth,
                                 'y_min_target':-pp.target_size_internal_grid_sigma*sigma_y_smooth,
                                 'y_max_target':pp.target_size_internal_grid_sigma*sigma_y_smooth,
                                 'Dh_target':.2*sigma_x_smooth},
                    save_pyecl_outp_as='cloud_evol_part%03d_ring%03d' % (self.SimSt.present_simulation_part, 
                        self.ring_of_CPUs.myring),
                    save_only=pp.save_only,
                    sparse_solver='PyKLU', enable_kick_x=pp.enable_kick_x, enable_kick_y=pp.enable_kick_y)
            print('Done.')



        # split the machine
        i_end_parallel = len(self.machine.one_turn_map)-self.n_non_parallelizable
        sharing = shs.ShareSegments(i_end_parallel, self.ring_of_CPUs.N_nodes_per_ring)
        i_start_part, i_end_part = sharing.my_part(self.ring_of_CPUs.myid_in_ring)
        self.mypart = self.machine.one_turn_map[i_start_part:i_end_part]

        if self.ring_of_CPUs.I_am_at_end_ring:
            self.non_parallel_part = self.machine.one_turn_map[i_end_parallel:]
            

        #install eclouds in my part
        if pp.enable_ecloud:
            my_new_part = []
            self.my_list_eclouds = []
            for ele in self.mypart:
                if ele in self.machine.transverse_map:
                    ecloud_new = ecloud.generate_twin_ecloud_with_shared_space_charge()
                    
                    # we save buildup info only for the first cloud in each ring
                    if self.ring_of_CPUs.myid_in_ring>0 or len(self.my_list_eclouds)>0:
                        ecloud_new.remove_savers()
                    
                    my_new_part.append(ecloud_new)
                    self.my_list_eclouds.append(ecloud_new)
                my_new_part.append(ele)

            self.mypart = my_new_part
            
            print(('Hello, I am %d.%d, my part looks like: %s. Saver status: %s'%(
                self.ring_of_CPUs.myring, self.ring_of_CPUs.myid_in_ring, self.mypart, 
                [(ec.cloudsim.cloud_list[0].pyeclsaver is not None) for ec in self.my_list_eclouds])))
            
       
    def init_master(self):
        
        import PyPARIS.gen_multibunch_beam as gmb
        from scipy.constants import c as clight, e as qe
        from PyHEADTAIL.particles.slicing import UniformBinSlicer
        
        if self.SimSt.first_run:
            if pp.load_beam_from_folder is None:
                print('Building the beam!')
                list_bunches = gmb.gen_matched_multibunch_beam(self.machine, pp.macroparticlenumber, 
                        pp.filling_pattern, pp.b_spac_s, pp.bunch_intensity, 
                        pp.epsn_x, pp.epsn_y, pp.sigma_z, pp.non_linear_long_matching, pp.min_inten_slice4EC)
                # compute and apply initial displacements
                inj_opt = self.machine.transverse_map.get_injection_optics()
                sigma_x = np.sqrt(inj_opt['beta_x']*pp.epsn_x/self.machine.betagamma)
                sigma_y = np.sqrt(inj_opt['beta_y']*pp.epsn_y/self.machine.betagamma)
                x_kick = pp.x_kick_in_sigmas*sigma_x
                y_kick = pp.y_kick_in_sigmas*sigma_y
                for bunch in list_bunches:
                    bunch.x += x_kick
                    bunch.y += y_kick
            else:
                # Load based on input
                list_bunches = gmb.load_multibunch_beam(pp.load_beam_from_folder)
        else:
            # Load from previous run
            print('Loading beam from file...')
            dirname = 'beam_status_part%02d'%(self.SimSt.present_simulation_part-1)
            list_bunches = gmb.load_multibunch_beam(dirname) 
            print('Loaded beam from file.')
        
        for bb in list_bunches:
            bb.slice_info['simstate_part'] = self.SimSt.present_simulation_part

        return list_bunches


    def init_start_ring(self):
        self.bunch_monitor = None
    
    def perform_bunch_operations_at_start_ring(self, bunch):
        
        if self.bunch_monitor is None:
            
            simstate_part = bunch.slice_info['simstate_part']
            stats_to_store = pp.stats_to_store

            n_stored_turns = np.sum(np.array(pp.filling_pattern)>0)*(\
                self.ring_of_CPUs.N_turns/self.ring_of_CPUs.N_parellel_rings\
                + self.ring_of_CPUs.N_parellel_rings)

            from PyHEADTAIL.monitors.monitors import BunchMonitor
            self.bunch_monitor = BunchMonitor(
                                'bunch_monitor_part%03d_ring%03d'%(
                                    simstate_part, self.ring_of_CPUs.myring),
                                n_stored_turns, 
                                {'Comment':'PyHDTL simulation'}, 
                                write_buffer_every = 1,
                                stats_to_store = stats_to_store)

            # define a slice monitor 
            z_left = bunch.slice_info['z_bin_left'] - bunch.slice_info['z_bin_center']
            z_right = bunch.slice_info['z_bin_right'] - bunch.slice_info['z_bin_center']
            slicer = UniformBinSlicer(n_slices = self.n_slices_per_bunch, z_cuts=(z_left, z_right))
            from PyHEADTAIL.monitors.monitors import SliceMonitor
            self.slice_monitor = SliceMonitor('slice_monitor_part%03d_ring%03d'%(
                simstate_part, self.ring_of_CPUs.myring),
                n_stored_turns, slicer,  {'Comment':'PyHDTL simulation'}, 
                write_buffer_every = 1, bunch_stats_to_store=stats_to_store,
                slice_stats_to_store=pp.slice_stats_to_store)

        # Save bunch properties
        if bunch.macroparticlenumber > 0 and bunch.slice_info['i_turn'] < self.N_turns:
            # Attach bound methods to monitor i_bunch and i_turns 
            bunch.i_bunch = types.MethodType(lambda ss: ss.slice_info['i_bunch'], bunch)
            bunch.i_turn = types.MethodType(lambda ss: ss.slice_info['i_turn'], bunch)
            self.bunch_monitor.dump(bunch)

            # Monitor slice wrt bunch center
            bunch.z -= bunch.slice_info['z_bin_center']
            self.slice_monitor.dump(bunch)
            bunch.z += bunch.slice_info['z_bin_center']

        # Save full beam at user-defined positions
        if bunch.slice_info['i_turn'] in self.save_beam_at_turns:
            dirname = 'bunch_states_turn%d'%bunch.slice_info['i_turn']
            import PyPARIS.gen_multibunch_beam as gmb
            gmb.save_bunch_to_folder(bunch, dirname)

        # Save full beam at end simulation 
        if bunch.slice_info['i_turn'] == self.N_turns:
            # PyPARIS wants N_turns to be a multiple of N_parellel_rings
            assert(self.ring_of_CPUs.I_am_the_master) 
            dirname = 'beam_status_part%02d'%(self.SimSt.present_simulation_part)
            import PyPARIS.gen_multibunch_beam as gmb
            gmb.save_bunch_to_folder(bunch, dirname)
            if bunch.slice_info['i_bunch'] == bunch.slice_info['N_bunches_tot_beam'] - 1:
                if not self.SimSt.first_run:
                    os.system('rm -r beam_status_part%02d' % (self.SimSt.present_simulation_part - 1))
                self.SimSt.after_simulation()

    def slice_bunch_at_start_ring(self, bunch):
        list_slices = sl.slice_a_bunch(bunch, self.z_cut_slicing, self.n_slices_per_bunch)
        return list_slices

    def treat_piece(self, piece):
        for ele in self.mypart: 
            ele.track(piece)
        
    def merge_slices_at_end_ring(self, list_slices):
        bunch = sl.merge_slices_into_bunch(list_slices)
        return bunch

    def perform_bunch_operations_at_end_ring(self, bunch):
        #finalize present turn (with non parallel part, e.g. synchrotron motion)
        if bunch.macroparticlenumber>0:
            for ele in self.non_parallel_part:
                ele.track(bunch)

    def piece_to_buffer(self, piece):
        buf = ch.beam_2_buffer(piece)
        return buf

    def buffer_to_piece(self, buf):
        piece = ch.buffer_2_beam(buf)
        return piece


class DummyComm(object):

    def __init__(self, N_cores_pretend, pretend_proc_id):
        self.N_cores_pretend = N_cores_pretend
        self.pretend_proc_id = pretend_proc_id

    def Get_size(self):
        return self.N_cores_pretend

    def Get_rank(self):
        return self.pretend_proc_id

    def Barrier(self):
        pass


def get_sim_instance(N_cores_pretend, id_pretend, 
                     init_sim_objects_auto=True):

    from PyPARIS.ring_of_CPUs_multiturn import RingOfCPUs_multiturn
    myCPUring = RingOfCPUs_multiturn(Simulation(),
            comm=DummyComm(N_cores_pretend, id_pretend))
    return myCPUring.sim_content



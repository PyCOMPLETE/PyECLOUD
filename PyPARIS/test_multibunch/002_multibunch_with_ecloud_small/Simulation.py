# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

import sys, os
BIN = os.path.expanduser("../../../")
sys.path.append(BIN)

import types

import numpy as np
from scipy.constants import c

import PyPARIS.communication_helpers as ch
import PyPARIS.share_segments as shs
import PyPARIS.slicing_tool as sl

verbose = False

sigma_z_bunch = 10e-2

machine_configuration = 'HLLHC-injection'
n_segments = 8

Qp_x = 0.
Qp_y = 0.

flag_aperture = True

enable_transverse_damper = False
dampingrate_x = 10000.
dampingrate_y = 50.

# Beam properties
non_linear_long_matching = False

bunch_intensity = 1e11
epsn_x = 2.5e-6
epsn_y = 2.5e-6
sigma_z = sigma_z_bunch

#Filling pattern: here head is left and tail is right
b_spac_s = 25e-9/5
filling_pattern = 4*([1.]+4*[0.])

macroparticlenumber = 1000000
min_inten_slice4EC = 1e7

x_kick_in_sigmas = 0.25
y_kick_in_sigmas = 0.25

target_size_internal_grid_sigma = 10.

enable_ecloud = True

enable_kick_x = True
enable_kick_y = False

L_ecloud_tot = 20e3

pickle_beam = False


class Simulation(object):
    def __init__(self):
        self.N_turns = 10000
        self.N_buffer_float_size = 10000000
        self.N_buffer_int_size = 20
        self.N_parellel_rings = 2
        
        self.n_slices_per_bunch = 200
        self.z_cut_slicing = 3*sigma_z_bunch
        self.N_pieces_per_transfer = 300
        

    def init_all(self):
        
        print('Exec init...')
        
        self.ring_of_CPUs.verbose = verbose

        from LHC_custom import LHC
        self.machine = LHC(n_segments = n_segments, machine_configuration = machine_configuration,
                        Qp_x=Qp_x, Qp_y=Qp_y)
        self.n_non_parallelizable = 1 #RF

        inj_optics = self.machine.transverse_map.get_injection_optics()
        sigma_x_smooth = np.sqrt(inj_optics['beta_x']*epsn_x/self.machine.betagamma)
        sigma_y_smooth = np.sqrt(inj_optics['beta_y']*epsn_y/self.machine.betagamma)

        if flag_aperture:
            # setup transverse losses (to "protect" the ecloud)
            import PyHEADTAIL.aperture.aperture as aperture
            apt_xy = aperture.EllipticalApertureXY(x_aper=target_size_internal_grid_sigma*sigma_x_smooth, 
                                                   y_aper=target_size_internal_grid_sigma*sigma_x_smooth)
            self.machine.one_turn_map.append(apt_xy)
            self.n_non_parallelizable +=1 

        if enable_transverse_damper:
            # setup transverse damper
            from PyHEADTAIL.feedback.transverse_damper import TransverseDamper
            damper = TransverseDamper(dampingrate_x=dampingrate_x, dampingrate_y=dampingrate_y)
            self.machine.one_turn_map.append(damper)
            self.n_non_parallelizable +=1
            
        if enable_ecloud:
            print('Build ecloud...')
            import PyECLOUD.PyEC4PyHT as PyEC4PyHT
            ecloud = PyEC4PyHT.Ecloud(
                    L_ecloud=L_ecloud_tot/n_segments, slicer=None, slice_by_slice_mode=True,
                    Dt_ref=5e-12, pyecl_input_folder='./pyecloud_config',
                    chamb_type = 'polyg' ,
                    filename_chm= 'LHC_chm_ver.mat', 
                    #init_unif_edens_flag=1,
                    #init_unif_edens=1e7,
                    #N_mp_max = 3000000,
                    #nel_mp_ref_0 = 1e7/(0.7*3000000),
                    #B_multip = [0.],
                    #~ PyPICmode = 'ShortleyWeller_WithTelescopicGrids',
                    #~ f_telescope = 0.3,
                    target_grid = {'x_min_target':-target_size_internal_grid_sigma*sigma_x_smooth, 'x_max_target':target_size_internal_grid_sigma*sigma_x_smooth,
                                   'y_min_target':-target_size_internal_grid_sigma*sigma_y_smooth,'y_max_target':target_size_internal_grid_sigma*sigma_y_smooth,
                                   'Dh_target':.2*sigma_x_smooth},
                    #~ N_nodes_discard = 10.,
                    #~ N_min_Dh_main = 10,
                    #x_beam_offset = x_beam_offset,
                    #y_beam_offset = y_beam_offset,
                    #probes_position = probes_position,
                    save_pyecl_outp_as = 'cloud_evol_ring%d'%self.ring_of_CPUs.myring,
                    sparse_solver = 'PyKLU', enable_kick_x=enable_kick_x, enable_kick_y=enable_kick_y)
            print('Done.')



        # split the machine
        i_end_parallel = len(self.machine.one_turn_map)-self.n_non_parallelizable
        sharing = shs.ShareSegments(i_end_parallel, self.ring_of_CPUs.N_nodes_per_ring)
        i_start_part, i_end_part = sharing.my_part(self.ring_of_CPUs.myid_in_ring)
        self.mypart = self.machine.one_turn_map[i_start_part:i_end_part]

        if self.ring_of_CPUs.I_am_at_end_ring:
            self.non_parallel_part = self.machine.one_turn_map[i_end_parallel:]
            

        #install eclouds in my part
        if enable_ecloud:
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
            
            print(('Hello, I am %d.%d, my part looks like: %s. Saver status: %s'%(self.ring_of_CPUs.myring, self.ring_of_CPUs.myid_in_ring, self.mypart, [(ec.cloudsim.cloud_list[0].pyeclsaver is not None) for ec in self.my_list_eclouds])))
            


        
    def init_master(self):
        
        print('Building the beam!')
        
        from scipy.constants import c as clight, e as qe
        from PyHEADTAIL.particles.slicing import UniformBinSlicer
        
        import PyPARIS.gen_multibunch_beam as gmb
        list_bunches = gmb.gen_matched_multibunch_beam(self.machine, macroparticlenumber, filling_pattern, b_spac_s, 
                bunch_intensity, epsn_x, epsn_y, sigma_z, non_linear_long_matching, min_inten_slice4EC)
                
        if pickle_beam:
            import pickle
            with open('init_beam.pkl', 'w') as fid:
                pickle.dump({'list_bunches': list_bunches}, fid)

        # compute and apply initial displacements
        inj_opt = self.machine.transverse_map.get_injection_optics()
        sigma_x = np.sqrt(inj_opt['beta_x']*epsn_x/self.machine.betagamma)
        sigma_y = np.sqrt(inj_opt['beta_y']*epsn_y/self.machine.betagamma)
        x_kick = x_kick_in_sigmas*sigma_x
        y_kick = y_kick_in_sigmas*sigma_y
        for bunch in list_bunches:
            bunch.x += x_kick
            bunch.y += y_kick


        return list_bunches

    def init_start_ring(self):
        stats_to_store = [
         'mean_x', 'mean_xp', 'mean_y', 'mean_yp', 'mean_z', 'mean_dp',
         'sigma_x', 'sigma_y', 'sigma_z','sigma_dp', 'epsn_x', 'epsn_y',
         'epsn_z', 'macroparticlenumber',
         'i_bunch', 'i_turn']

        n_stored_turns = len(filling_pattern)*(self.ring_of_CPUs.N_turns/self.ring_of_CPUs.N_parellel_rings + self.ring_of_CPUs.N_parellel_rings)

        from PyHEADTAIL.monitors.monitors import BunchMonitor
        self.bunch_monitor = BunchMonitor('bunch_monitor_ring%03d'%self.ring_of_CPUs.myring,
                            n_stored_turns, 
                            {'Comment':'PyHDTL simulation'}, 
                            write_buffer_every = 1,
                            stats_to_store = stats_to_store)

    def perform_bunch_operations_at_start_ring(self, bunch):
        # Attach bound methods to monitor i_bunch and i_turns 
        # (In the future we might upgrade PyHEADTAIL to pass the lambda to the monitor)
        if bunch.macroparticlenumber>0:
            bunch.i_bunch = types.MethodType(lambda self: self.slice_info['i_bunch'], bunch)
            bunch.i_turn = types.MethodType(lambda self: self.slice_info['i_turn'], bunch)
            self.bunch_monitor.dump(bunch)        

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





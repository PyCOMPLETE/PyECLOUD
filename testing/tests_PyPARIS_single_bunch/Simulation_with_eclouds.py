import sys, os
BIN = os.path.expanduser("../../")
sys.path.append(BIN)

import numpy as np
from scipy.constants import c, e

import PyPARIS.share_segments as shs
import PyPARIS.communication_helpers as ch

# compute sigma x and y
epsn_x = 2.5e-6
epsn_y = 2.5e-6
B_multip = [0.5]
n_slices = 64
n_segments = 3

filename = '../tests_PyEC4PyHT/headtail_for_test/test_protons/SPS_Q20_proton_check_dipole_3kicks_20150212_prb.dat'
B_multip = [0.5]
N_turns_to_simulate = 8


class Simulation(object):
	def __init__(self):
		self.N_turns = N_turns_to_simulate
		self.N_pieces_per_transfer=10

	def init_all(self):

		
		self.n_slices = n_slices
		self.n_segments = n_segments

		from machines_for_testing import SPS
		self.machine = SPS(n_segments = n_segments, 
			machine_configuration = 'Q20-injection', accQ_x=20., accQ_y=20., 
					RF_at='end_of_transverse')

		
		# We suppose that all the object that cannot be slice parallelized are at the end of the ring
		i_end_parallel = len(self.machine.one_turn_map)-1 #only RF is not parallelizable

		# split the machine
		sharing = shs.ShareSegments(i_end_parallel, self.ring_of_CPUs.N_nodes)
		myid = self.ring_of_CPUs.myid
		i_start_part, i_end_part = sharing.my_part(myid)
		self.mypart = self.machine.one_turn_map[i_start_part:i_end_part]
		if self.ring_of_CPUs.I_am_a_worker:
			print 'I am id=%d (worker) and my part is %d long'%(myid, len(self.mypart))
		elif self.ring_of_CPUs.I_am_the_master:
			self.non_parallel_part = self.machine.one_turn_map[i_end_parallel:]
			print 'I am id=%d (master) and my part is %d long'%(myid, len(self.mypart))

	
		# config e-cloud
		init_unif_edens_flag=1
		init_unif_edens=2e11
		N_MP_ele_init = 100000
		N_mp_max = N_MP_ele_init*4.
		
		# define apertures and Dh_sc to simulate headtail 
		inj_optics = self.machine.transverse_map.get_injection_optics()
		sigma_x = np.sqrt(inj_optics['beta_x']*epsn_x/self.machine.betagamma)
		sigma_y = np.sqrt(inj_optics['beta_y']*epsn_y/self.machine.betagamma)
		x_aper  = 20*sigma_x
		y_aper  = 20*sigma_y
		Dh_sc = 2*x_aper/128/2
		
		# initial MP size
		nel_mp_ref_0 = init_unif_edens*4*x_aper*y_aper/N_MP_ele_init

		import PyECLOUD.PyEC4PyHT as PyEC4PyHT
		ecloud = PyEC4PyHT.Ecloud(slice_by_slice_mode=True,
						L_ecloud=self.machine.circumference/n_segments, 
						slicer=None, 
						Dt_ref=25e-12, 
						pyecl_input_folder='../tests_PyEC4PyHT/drift_sim/',
						x_aper=x_aper, y_aper=y_aper, Dh_sc=Dh_sc,
						init_unif_edens_flag=init_unif_edens_flag,
						init_unif_edens=init_unif_edens, 
						N_mp_max=N_mp_max,
						nel_mp_ref_0=nel_mp_ref_0,
						B_multip=B_multip)
		
		
		my_new_part = []
		self.my_list_eclouds = []
		for ele in self.mypart:
			my_new_part.append(ele)
			if ele in self.machine.transverse_map:
				ecloud_new = ecloud.generate_twin_ecloud_with_shared_space_charge()
				my_new_part.append(ecloud_new)
				self.my_list_eclouds.append(ecloud_new)
		self.mypart = my_new_part

	def init_master(self):
		
		# beam parameters
		sigma_z = 0.2
		intensity = 1.15e11
		macroparticlenumber_track = 300000

		# initialization bunch
		bunch = self.machine.generate_6D_Gaussian_bunch(
			macroparticlenumber_track, intensity, epsn_x, epsn_y, sigma_z=sigma_z)
		print 'Bunch initialized.'

		#replace particles with HDTL ones
		self.n_part_per_turn = 5000
		appo = np.loadtxt(filename)
		
		parid = np.reshape(appo[:,0], (-1, self.n_part_per_turn))[::self.n_segments,:]
		x = np.reshape(appo[:,1], (-1, self.n_part_per_turn))[::self.n_segments,:]
		xp = np.reshape(appo[:,2], (-1, self.n_part_per_turn))[::self.n_segments,:]
		y = np.reshape(appo[:,3], (-1, self.n_part_per_turn))[::self.n_segments,:]
		yp =np.reshape(appo[:,4], (-1, self.n_part_per_turn))[::self.n_segments,:]
		z = np.reshape(appo[:,5], (-1, self.n_part_per_turn))[::self.n_segments,:]
		zp = np.reshape(appo[:,6], (-1, self.n_part_per_turn))[::self.n_segments,:]

		# replace first particles with HEADTAIL ones
		bunch.x[:self.n_part_per_turn] = x[0,:]
		bunch.xp[:self.n_part_per_turn] = xp[0,:]
		bunch.y[:self.n_part_per_turn] = y[0,:]
		bunch.yp[:self.n_part_per_turn] = yp[0,:]
		bunch.z[:self.n_part_per_turn] = z[0,:]
		bunch.dp[:self.n_part_per_turn] =zp[0,:]

		# save id and momenta before track
		self.id_before = bunch.id[bunch.id<=self.n_part_per_turn]
		self.xp_before = bunch.xp[bunch.id<=self.n_part_per_turn]
		self.yp_before = bunch.yp[bunch.id<=self.n_part_per_turn]

		# initial slicing
		from PyHEADTAIL.particles.slicing import UniformBinSlicer
		self.slicer = UniformBinSlicer(n_slices = self.n_slices, n_sigma_z = 3.)

		self.rms_err_x_list = []
		self.rms_err_y_list = []
		
		#slice for the first turn
		slice_obj_list = bunch.extract_slices(self.slicer)

		pieces_to_be_treated = slice_obj_list
		
		print 'N_turns', self.N_turns

		return pieces_to_be_treated

	def init_worker(self):
		pass

	def treat_piece(self, piece):
		for ele in self.mypart: 
				ele.track(piece)

	def finalize_turn_on_master(self, pieces_treated):
		
		# re-merge bunch
		bunch = sum(pieces_treated)

		#finalize present turn (with non parallel part, e.g. synchrotron motion)
		for ele in self.non_parallel_part:
			ele.track(bunch)

		# id and momenta after track
		id_after = bunch.id[bunch.id<=self.n_part_per_turn]
		xp_after = bunch.xp[bunch.id<=self.n_part_per_turn]
		z_after = bunch.z[bunch.id<=self.n_part_per_turn]
		yp_after = bunch.yp[bunch.id<=self.n_part_per_turn]

		# sort id and momenta after track
		indsort = np.argsort(id_after)
		id_after = np.take(id_after, indsort)
		xp_after = np.take(xp_after, indsort)
		yp_after = np.take(yp_after, indsort)
		z_after = np.take(z_after, indsort)

		# save results
		import PyPARIS.myfilemanager as mfm
		mfm.dict_to_h5({\
		    'id_after': id_after,
		    'xp_after': xp_after,
		    'yp_after': yp_after,
		    'z_after': z_after,
		    'id_before':self.id_before,
		    'xp_before':self.xp_before,
		    'yp_before':self.yp_before},
                    'particles_at_turn_%d.h5'%self.ring_of_CPUs.i_turn)


		# prepare next turn (re-slice)
		new_pieces_to_be_treated = bunch.extract_slices(self.slicer)
		orders_to_pass = ['reset_clouds']

		return orders_to_pass, new_pieces_to_be_treated


	def execute_orders_from_master(self, orders_from_master):
		if 'reset_clouds' in orders_from_master:
			for ec in self.my_list_eclouds: ec.finalize_and_reinitialize()


		
	def finalize_simulation(self):
		pass
	def piece_to_buffer(self, piece):
		buf = ch.beam_2_buffer(piece)
		return buf

	def buffer_to_piece(self, buf):
		piece = ch.buffer_2_beam(buf)
		return piece





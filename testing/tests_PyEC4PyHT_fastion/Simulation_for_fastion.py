import sys
sys.path.append("../../../")

import communication_helpers as ch
import share_segments as shs

import time
import numpy as np
import scipy.io as sio
import pickle
from scipy.constants import c, e, m_p


print 'Start initialization'
t_start = time.mktime(time.localtime())

chamb_type = 'rect'
x_aper = 5.e-4
y_aper = 5.e-5
filename_chm = None
B_multip_per_eV = [0.]
B_multip_per_eV = np.array(B_multip_per_eV)

# define PIC grid size
Dh_x = x_aper/125.
Dh_y = y_aper/125.
Dh_sc = [Dh_x, Dh_y]

init_unif_edens_flag = -1
init_unif_edens = None

gas_ion_flag = 1
unif_frac = 0. 
P_nTorr = 20. 
sigma_ion_MBarn = 1.5
Temp_K = 300.
A = 18
E_init_ion = 0.157e-27 * A
ion_mass = A * m_p
ion_charge = e

#N_MP_ele_init = 501*312
N_MP_ele_init = 501 * 200
#N_MP_ele_init = 501*156
#N_MP_ele_init = 501*312 / 2.
N_mp_max = N_MP_ele_init

machine_configuration ='CLIC_DR_2GHz'

# define bunch
intensity = 4.1e9
n_macroparticles = 4688
#n_macroparticles = 10000 * 10

if machine_configuration =='CLIC_DR_1GHz':
	sigma_z = 1.8e-3
	epsn_x = 0.456e-6
elif machine_configuration =='CLIC_DR_2GHz':
	sigma_z = 1.6e-3
	epsn_x = 0.472e-6
	epsn_x = 0.01e-6
epsn_y = 0.0048e-6

optics_mode = 'smooth'
#optics_mode = 'non-smooth'
#optics_file = 'CLIC_DR_n260_optics.pkl'
#optics_file = 'CLIC_DR_n677_dist_55_cm_optics.pkl'

#n_bunches = 312
n_bunches = 64
#n_segments = 260.
n_segments = 79.
#n_segments = 677. 
#n_segments = 339.
N_turns = 1



class Simulation(object):
	def __init__(self):
		self.N_turns = N_turns
		self.N_pieces_per_transfer = 1
		self.N_buffer_float_size = 75000

	def init_all(self):

		self.n_slices = n_bunches
		self.n_segments = n_segments

		# define the machine
		from CLIC_DR import CLIC_DR
		if optics_mode == 'smooth':
			self.machine = CLIC_DR(machine_configuration=machine_configuration, n_segments=n_segments)

		elif optics_mode == 'non-smooth':
			with open(optics_file) as fid:
				optics = pickle.load(fid)
			optics.pop('circumference')
			# optics_new = {}
			# for kk in optics.keys():
			# 	optics_new[kk] = optics[kk][::2]
			# 	np.append(optics_new[kk], optics[kk][-1])
			# optics = optics_new
			self.machine = CLIC_DR(machine_configuration=machine_configuration, optics_mode = 'non-smooth',  **optics)		

		# define MP size
		#nel_mp_ref_0 = P_nTorr * sigma_ion_MBarn / 37.89
		#nel_mp_ref_0 = P_nTorr * sigma_ion_MBarn / 37.89 * 2
		#nel_mp_ref_0 = P_nTorr * sigma_ion_MBarn / 37.89 / 156.
		nel_mp_ref_0 = P_nTorr * sigma_ion_MBarn / 37.89 / 200.

		# prepare e-cloud
		import PyECLOUD.PyEC4PyHT_fastion as PyEC4PyHT_fi
		ecloud = PyEC4PyHT_fi.Ecloud_fastion(L_ecloud = self.machine.circumference/n_segments, slicer = None, 
												Dt_ref = 1e-9, MP_e_mass = ion_mass, MP_e_charge = ion_charge, 
												pyecl_input_folder = './pyecloud_config', 
												slice_by_slice_mode = True, 
												#beam_monitor = beam_monitor,
												ionize_only_first_bunch = True,
												chamb_type = chamb_type, PyPICmode = 'FFT_OpenBoundary',
												x_aper = x_aper, y_aper = y_aper,
												filename_chm = filename_chm, Dh_sc = Dh_sc,
												init_unif_edens_flag = init_unif_edens_flag,
												init_unif_edens = init_unif_edens, 
												gas_ion_flag = gas_ion_flag, unif_frac = unif_frac, 
												P_nTorr = P_nTorr, sigma_ion_MBarn = sigma_ion_MBarn, 
												Temp_K = Temp_K, E_init_ion = E_init_ion,
												N_mp_max = N_mp_max,
												nel_mp_ref_0 = nel_mp_ref_0,
												B_multip = B_multip_per_eV*self.machine.p0/e*c,
												switch_model = 'perfect_absorber')


		nx, ny = ecloud.spacech_ele.PyPICobj.nx, ecloud.spacech_ele.PyPICobj.ny
		print 'nx = %d, ny = %d'%(nx, ny)

		n_non_parallelizable = 2 #rf and wrapper
		
		# We suppose that all the object that cannot be slice parallelized are at the end of the ring
		i_end_parallel = len(self.machine.one_turn_map) - n_non_parallelizable

		# split the machine
		sharing = shs.ShareSegments(i_end_parallel, self.ring_of_CPUs.N_nodes)
		myid = self.ring_of_CPUs.myid
		i_start_part, i_end_part = sharing.my_part(myid)
		self.mypart = self.machine.one_turn_map[i_start_part:i_end_part]
		if self.ring_of_CPUs.I_am_a_worker:
			print 'I am id=%d/%d (worker) and my part is %d long'%(myid, self.ring_of_CPUs.N_nodes, len(self.mypart))
		elif self.ring_of_CPUs.I_am_the_master:
			self.non_parallel_part = self.machine.one_turn_map[i_end_parallel:]
			print 'I am id=%d/%d (master) and my part is %d long'%(myid, self.ring_of_CPUs.N_nodes, len(self.mypart))
		
		#install eclouds in my part
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

		# beam position from fastion
		# fastion = np.loadtxt('test_hdtl2.dat')
		# x_fi = fastion[:,0]
		# y_fi = fastion[:,1]
		# xp_fi = fastion[:,2]
		# yp_fi = fastion[:,3]
		
		# generate beam
		print 'Initializing', n_bunches, 'bunches, of', n_macroparticles, 'macroparticles'
		print 'Using initial beam conditions from FASTION'

		bunches = []
		for i_bun in xrange(n_bunches):
			print 'Bunch', i_bun
			bunch = self.machine.generate_6D_Gaussian_bunch(n_macroparticles=n_macroparticles, intensity=intensity, epsn_x=epsn_x, epsn_y=epsn_y, sigma_z=sigma_z)
			bunch.z -= self.machine.circumference / self.machine.longitudinal_map.harmonics[0]*i_bun

			
			# bunch.x = x_fi[1e4*i_bun:1e4*(i_bun+1)]
			# bunch.y = y_fi[1e4*i_bun:1e4*(i_bun+1)]
			# bunch.xp = xp_fi[1e4*i_bun:1e4*(i_bun+1)]
			# bunch.yp = yp_fi[1e4*i_bun:1e4*(i_bun+1)]

			mean_x = bunch.mean_x()
			mean_y = bunch.mean_y()	

			print 'Bunch centroid at', bunch.mean_x(), bunch.mean_y()
			bunches.append(bunch)
	
		beam = sum(bunches)


		# initial slicing
		self.slicer = self.machine.buncher
		slices = beam.get_slices(self.slicer)
		self.slicer.add_statistics(sliceset=slices, beam=beam, statistics=True)

		mask_filled = slices.n_macroparticles_per_slice > 0
		n_filled = np.sum(mask_filled)
		z_cuts_filled = (np.min(slices.z_bins[mask_filled]), np.max(slices.z_bins[np.append(mask_filled, True)]))

		from PyHEADTAIL.particles.slicing import UniformBinSlicer
		filled_buncher = UniformBinSlicer(n_filled, z_cuts=z_cuts_filled)
		filled_slices = beam.get_slices(filled_buncher)
		filled_buncher.add_statistics(sliceset=filled_slices, beam=beam, statistics=True)

		# define a beam monitor 
		from PyHEADTAIL.monitors.monitors import SliceMonitor
		self.beam_monitor = SliceMonitor(filename='PyEC4PyHT_parallel_A%d_%db_%dips_%dturns_%.2fnTorr'%(A, n_bunches, n_segments, N_turns, P_nTorr), 
											n_steps = N_turns, slicer=filled_buncher, write_buffer_every=100)


		#slice for the first turn
		slice_obj_list = beam.extract_slices(self.slicer)

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
		beam = sum(pieces_treated)

		#finalize present turn (with non parallel part, e.g. synchrotron motion)
		for ele in self.non_parallel_part:
			ele.track(beam)
			
		# save results		
		#print '%s Turn %d'%(time.strftime("%d/%m/%Y %H:%M:%S", time.localtime()), i_turn)
		self.beam_monitor.dump(beam)
		
		# prepare next turn (re-slice)
		new_pieces_to_be_treated = beam.extract_slices(self.slicer)
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


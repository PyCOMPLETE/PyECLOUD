import sys, os
sys.path.append(os.path.expanduser('../../../'))
sys.path.append(os.path.expanduser('../../../PyHEADTAIL/'))


import numpy as np
import pylab as pl
import PyECLOUD.mystyle as ms
import time


B_multip = [0.5]
N_kicks = 3
N_turns = 2


pl.close('all')
ms.mystyle(fontsz=14)


# define machine for PyHEADTAIL
from PyHEADTAIL.particles.slicing import UniformBinSlicer
from machines_for_testing import SPS
machine = SPS(n_segments=N_kicks, machine_configuration='Q20-injection', accQ_x=20., accQ_y=20.)
machine.one_turn_map.remove(machine.longitudinal_map) # We apply it separately


# compute sigma x and y
epsn_x = 2.5e-6
epsn_y = 2.5e-6

inj_optics = machine.transverse_map.get_injection_optics()
sigma_x = np.sqrt(inj_optics['beta_x'] * epsn_x / machine.betagamma)
sigma_y = np.sqrt(inj_optics['beta_y'] * epsn_y / machine.betagamma)

# define apertures and Dh_sc to simulate headtail conditions
x_aper = 20 * sigma_x
y_aper = 20 * sigma_y
Dh_sc = 2 * x_aper / 128 / 2

# ecloud
import PyECLOUD.PyEC4PyHT as PyEC4PyHT
from PyHEADTAIL.particles.slicing import UniformBinSlicer
slicer = UniformBinSlicer(n_slices=64, n_sigma_z=3.)

init_unif_edens_flag = 1
init_unif_edens = 2e11
N_MP_ele_init = 100000
N_mp_max = N_MP_ele_init * 4.

nel_mp_ref_0 = init_unif_edens * 4 * x_aper * y_aper / N_MP_ele_init

new_one_turn_map = []
ecloud_list = []
for ele in machine.one_turn_map:
	new_one_turn_map.append(ele)
	if ele in machine.transverse_map:
		new_ecloud = PyEC4PyHT.Ecloud(L_ecloud=machine.circumference / N_kicks, slicer=slicer,
                                Dt_ref=25e-12, pyecl_input_folder='./drift_sim',
                                x_aper=x_aper, y_aper=y_aper, Dh_sc=Dh_sc,
                                init_unif_edens_flag=init_unif_edens_flag,
                                init_unif_edens=init_unif_edens,
                                N_mp_max=N_mp_max,
                                nel_mp_ref_0=nel_mp_ref_0,
                                B_multip=B_multip, slice_by_slice_mode=True)
		new_one_turn_map.append(new_ecloud)
		ecloud_list.append(new_ecloud)

machine.one_turn_map = new_one_turn_map

# generate a bunch
bunch = machine.generate_6D_Gaussian_bunch(n_macroparticles=30000, intensity=1.15e11, epsn_x=epsn_x, epsn_y=epsn_y, sigma_z=0.2)


t_start_slice = time.mktime(time.localtime())
for ii in range(N_turns):
	slices_list = bunch.extract_slices(slicer)

	for slice_obj in slices_list[::-1]:
		machine.track(slice_obj)  # , verbose = True)
	print('Turn', ii)

	bunch = sum(slices_list)

	machine.longitudinal_map.track(bunch)

	for ec in ecloud_list:
		ec.finalize_and_reinitialize()
t_end_slice = time.mktime(time.localtime())
print('Sliced %.2e s per turn'%((t_end_slice - t_start_slice) / float(N_turns)))

# Simulate bunch mode
machine_whole_bunch = SPS(n_segments=N_kicks, machine_configuration='Q20-injection', accQ_x=20., accQ_y=20.)


ecloud = PyEC4PyHT.Ecloud(L_ecloud=machine.circumference / N_kicks, slicer=slicer,
                          Dt_ref=25e-12, pyecl_input_folder='./drift_sim',
                          x_aper=x_aper, y_aper=y_aper, Dh_sc=Dh_sc,
                          init_unif_edens_flag=init_unif_edens_flag,
                          init_unif_edens=init_unif_edens,
                          N_mp_max=N_mp_max,
                          nel_mp_ref_0=nel_mp_ref_0,
                          B_multip=B_multip, slice_by_slice_mode=False)

machine_whole_bunch.install_after_each_transverse_segment(ecloud)

t_start_bunch = time.mktime(time.localtime())
for ii in range(N_turns):
	print('Turn', ii)
	machine_whole_bunch.track(bunch)
t_end_bunch = time.mktime(time.localtime())

print('\n\n')
print('Sliced %.2e s per turn'%((t_end_slice - t_start_slice) / float(N_turns)))
print('Full bunch %.2e s per turn'%((t_end_bunch - t_start_bunch) / float(N_turns)))

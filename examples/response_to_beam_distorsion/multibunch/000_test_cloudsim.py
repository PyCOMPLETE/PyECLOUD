import sys
import os
BIN = os.path.expanduser("../../../../")  # folder containing PyECLOUD, PyPIC, PyKLU, PyHEADTAIL, PyPARIS
if BIN not in sys.path:
    sys.path.append(BIN)

from scipy.constants import c as clight
import numpy as np

import PyPARIS.gen_multibunch_beam as gmb
from PyHEADTAIL.particles.slicing import UniformBinSlicer

from machines_for_testing import LHC


n_segments = 1
bunch_intensity = 1e11
epsn_x = 2.5e-6
epsn_y = 2.5e-6
sigma_z = 1.000000e-09 / 4. * 299792458.

machine_configuration = 'Injection'

n_slices = 150
z_cut = 2.5e-9 * clight / 2

min_inten_slice4EC = 1e3
non_linear_long_matching = False
b_spac_s = 25e-9
#Here head is left and tail is right
filling_pattern = 1 * (60 * [1.] + 5 * [0])

macroparticlenumber = 100000

machine = LHC(machine_configuration=machine_configuration, beta_x=85.00, beta_y=90.0,
              optics_mode='smooth', n_segments=n_segments, RF_at='end_of_transverse')

list_bunches = gmb.gen_matched_multibunch_beam(machine, macroparticlenumber, filling_pattern, b_spac_s, bunch_intensity, epsn_x, epsn_y, sigma_z, non_linear_long_matching, min_inten_slice4EC)

for bb in list_bunches[::-1][30:]:
    bb.x += 3e-3

# Slice bunches
import PyPARIS.slicing_tool as st
list_slices = []
for bb in list_bunches:
    these_slices = st.slice_a_bunch(bb, z_cut=z_cut, n_slices=n_slices)
    list_slices += these_slices

# Build e-cloud
print('Build ecloud...')
import PyECLOUD.PyEC4PyHT as PyEC4PyHT
ecloud = PyEC4PyHT.Ecloud(
    L_ecloud=1., slicer=None, slice_by_slice_mode=True,
    Dt_ref=5e-12, pyecl_input_folder='./pyecloud_config',
    chamb_type='polyg' ,
    filename_chm='LHC_chm_ver.mat',
    #init_unif_edens_flag=1,
        #init_unif_edens=1e7,
        #N_mp_max = 3000000,
        #nel_mp_ref_0 = 1e7/(0.7*3000000),
        #B_multip = [0.],
        #~ PyPICmode = 'ShortleyWeller_WithTelescopicGrids',
        #~ f_telescope = 0.3,
        target_grid={'x_min_target': -5 * list_bunches[-1].sigma_x(), 'x_max_target': 5 * list_bunches[-1].sigma_x(),
                     'y_min_target': -5 * list_bunches[-1].sigma_y(), 'y_max_target': 5 * list_bunches[-1].sigma_y(),
                     'Dh_target': .2 * list_bunches[-1].sigma_x()},
        #~ N_nodes_discard = 10.,
        #~ N_min_Dh_main = 10,
        #x_beam_offset = x_beam_offset,
        #y_beam_offset = y_beam_offset,
        #probes_position = probes_position,
        save_pyecl_outp_as='test_saving',
        sparse_solver='PyKLU')
print('Done.')


# REMEBMBER TO START POPPING FROM THE RIGHT SIDE
print('Start cloud sim')
for ii in range(len(list_slices) - 1, -1, -1):
    ecloud.track(list_slices[ii])
#~ for cc in ecloud.cloudsim.cloud_list:
    #~ cc.pyeclsaver.t_last_save = 0.
    #~ cc.MP_e.nel_mp_ref = cc.MP_e.nel_mp_ref_0

print('End cloud sim')

# Some plotting
bucket_length_m = machine.circumference / (machine.longitudinal_map.harmonics[0])
b_spac_m = b_spac_s * machine.beta * clight
b_spac_buckets = np.round(b_spac_m / bucket_length_m)

beam = sum(list_bunches)

# Build profile of the full beam
thin_slicer = UniformBinSlicer(n_slices=10000, z_cuts=(-len(filling_pattern) * bucket_length_m * b_spac_buckets, bucket_length_m))
thin_slice_set = beam.get_slices(thin_slicer, statistics=True)

import matplotlib.pyplot as plt

plt.close('all')
plt.figure(1)


#~ import json
sp1 = plt.subplot(3, 1, 1)
sp2 = plt.subplot(3, 1, 2, sharex=sp1)
sp3 = plt.subplot(3, 1, 3, sharex=sp1)
sp1.plot(thin_slice_set.z_centers, thin_slice_set.charge_per_slice)

for ibuf, ss in enumerate(list_slices):
        sp1.axvline(x=ss.slice_info['z_bin_center'], color='k', alpha=0.5, linestyle='--')
        sp1.axvspan(xmin=ss.slice_info['z_bin_left'], xmax=ss.slice_info['z_bin_right'],
                    color={0: 'r', 1: 'b'}[ibuf%2], alpha=0.3)
        sp2.stem([ss.slice_info['z_bin_center']], [ss.slice_info['interact_with_EC']])
        sp3.stem([ss.slice_info['z_bin_center']], [ss.slice_info['i_slice']])
sp2.grid('on')


plt.figure(2)
spb1 = plt.subplot(3, 1, 1, sharex=sp1)
spb1.plot(thin_slice_set.z_centers, thin_slice_set.charge_per_slice)
spb2 = plt.subplot(3, 1, 2, sharex=sp1)
spb3 = plt.subplot(3, 1, 3, sharex=sp1)
for ibun, bun in enumerate(list_bunches):
        spb1.axvline(x=bun.slice_info['z_bin_center'], color='k', alpha=0.5, linestyle='--')
        spb1.axvspan(xmin=bun.slice_info['z_bin_left'], xmax=bun.slice_info['z_bin_right'],
                     color={0: 'r', 1: 'b'}[ibun%2], alpha=0.3)
        spb2.stem([bun.slice_info['z_bin_center']], [bun.slice_info['interact_with_EC']])
        spb3.stem([bun.slice_info['z_bin_center']], [bun.slice_info['i_bunch']])
spb2.grid('on')


plt.show()

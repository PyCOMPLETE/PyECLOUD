import sys
import os
BIN = os.path.expanduser("../../../") #folder containing PyECLOUD, PyPIC, PyKLU, PyHEADTAIL, PyPARIS
if BIN not in sys.path:
    sys.path.append(BIN)

from scipy.constants import c as clight
import numpy as np

import PyPARIS.gen_multibunch_beam as gmb
from PyHEADTAIL.particles.slicing import UniformBinSlicer

from machines_for_testing import SPS


n_segments = 1
bunch_intensity = 1e11
epsn_x = 2.5e-6
epsn_y = 3.5e-6
sigma_z = 10e-2

n_slices = 20
z_cut = 2*sigma_z

min_inten_slice4EC = 1e3
non_linear_long_matching = False
b_spac_s = 25e-9
#Here head is left and tail is right
filling_pattern = 6*[1]

macroparticlenumber = 100000

machine = SPS(n_segments = n_segments, 
            machine_configuration = 'Q20-injection', accQ_x=20., accQ_y=20., 
            RF_at='end_of_transverse', longitudinal_mode = 'non-linear')
            
list_bunches = gmb.gen_matched_multibunch_beam(machine, macroparticlenumber, filling_pattern, b_spac_s, bunch_intensity, epsn_x, epsn_y, sigma_z, non_linear_long_matching, min_inten_slice4EC)

# Slice bunches
import PyPARIS.slicing_tool as st
list_slices = []
for bb in list_bunches:
    these_slices = st.slice_a_bunch(bb, z_cut=z_cut, n_slices=n_slices)
    list_slices+=these_slices

# REMEBMBER TO START POPPING FROM THE RIGHT SIDE



# Some plotting
bucket_length_m = machine.circumference/(machine.longitudinal_map.harmonics[0])
b_spac_m =  b_spac_s*machine.beta*clight
b_spac_buckets = np.round(b_spac_m/bucket_length_m)

beam = sum(list_bunches)

# Build profile of the full beam
thin_slicer = UniformBinSlicer(n_slices=10000, z_cuts=(-len(filling_pattern)*bucket_length_m*b_spac_buckets, bucket_length_m))
thin_slice_set = beam.get_slices(thin_slicer, statistics=True)

import matplotlib.pyplot as plt

plt.close('all')
plt.figure(1)


#~ import json
sp1 = plt.subplot(3,1,1)
sp2 = plt.subplot(3,1,2, sharex=sp1)
sp3 = plt.subplot(3,1,3, sharex=sp1)
sp1.plot(thin_slice_set.z_centers, thin_slice_set.charge_per_slice)

for ibuf, ss in enumerate(list_slices):
        sp1.axvline(x=ss.slice_info['z_bin_center'], color='k', alpha=0.5, linestyle='--')
        sp1.axvspan(xmin=ss.slice_info['z_bin_left'], xmax=ss.slice_info['z_bin_right'],
            color={0:'r', 1:'b'}[ibuf%2], alpha = 0.3)
        sp2.stem([ss.slice_info['z_bin_center']], [ss.slice_info['interact_with_EC']])
        sp3.stem([ss.slice_info['z_bin_center']], [ss.slice_info['i_slice']])
sp2.grid('on')


    
plt.figure(2)
spb1 = plt.subplot(3,1,1, sharex=sp1)
spb1.plot(thin_slice_set.z_centers, thin_slice_set.charge_per_slice)
spb2 = plt.subplot(3,1,2, sharex=sp1)
spb3 = plt.subplot(3,1,3, sharex=sp1)
for ibun, bun in enumerate(list_bunches):
        spb1.axvline(x=bun.slice_info['z_bin_center'], color='k', alpha=0.5, linestyle='--')
        spb1.axvspan(xmin=bun.slice_info['z_bin_left'], xmax=bun.slice_info['z_bin_right'],
            color={0:'r', 1:'b'}[ibun%2], alpha = 0.3)
        spb2.stem([bun.slice_info['z_bin_center']], [bun.slice_info['interact_with_EC']])
        spb3.stem([bun.slice_info['z_bin_center']], [bun.slice_info['i_bunch']])
spb2.grid('on')



plt.show()

# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

import sys
sys.path.append('../../../')

import numpy as np
from scipy.constants import c as clight, e as qe

from PyHEADTAIL.particles.slicing import UniformBinSlicer
import PyPARIS.gen_multibunch_beam as gmb
import PyECLOUD.mystyle as ms

from machines_for_testing import SPS

n_segments = 1
bunch_intensity = 1e11
epsn_x = 2.5e-6
epsn_y = 3.5e-6
sigma_z = 10e-2

n_slices = 20
z_cut = 3*sigma_z

min_inten_slice4EC = 1e3

non_linear_long_matching = False


# ~ #Here head is left and tail is right
b_spac_s = 25e-9
filling_pattern = [1., 1., 0., 1., 1.]

# ~ b_spac_s = 5e-9
# ~ filling_pattern = 5*([1.]+4*[0.])

macroparticlenumber = 500000

# Build machine
machine = SPS(n_segments = n_segments, 
            machine_configuration = 'Q20-injection', accQ_x=20., accQ_y=20., 
            RF_at='end_of_transverse', longitudinal_mode = 'non-linear')
            
bucket_length_m = machine.circumference/(machine.longitudinal_map.harmonics[0])
b_spac_m =  b_spac_s*machine.beta*clight
b_spac_buckets = np.round(b_spac_m/bucket_length_m)
            
#generate bunch objects
list_bunches = gmb.gen_matched_multibunch_beam(machine, macroparticlenumber, filling_pattern, b_spac_s, bunch_intensity, epsn_x, epsn_y, sigma_z, non_linear_long_matching, min_inten_slice4EC)

beam = sum(list_bunches)

import PyPARIS.slicing_tool as st
import PyPARIS.communication_helpers as ch

# Turn slices into buffer
list_buffers = []
list_bunch_buffers = []
for bb in list_bunches:
    these_slices = st.slice_a_bunch(bb, z_cut=z_cut, n_slices=n_slices)
    list_bunch_buffers.append([])
    for ss in these_slices:
        thisbuffer = ch.beam_2_buffer(ss,verbose=True, mode='pickle')
        list_buffers.append(thisbuffer)
        list_bunch_buffers[-1].append(thisbuffer)

big_buffer = ch.combine_float_buffers(list_buffers)


# Build profile of the full beam
thin_slicer = UniformBinSlicer(n_slices=10000, z_cuts=(-len(filling_pattern)*bucket_length_m*b_spac_buckets, bucket_length_m))
thin_slice_set = beam.get_slices(thin_slicer, statistics=True)

import matplotlib.pyplot as plt
plt.close('all')

ms.mystyle_arial(fontsz=14, dist_tick_lab=5)



# re-split buffer
list_buffers_rec = ch.split_float_buffers(big_buffer)

# Plot including sub-slicing
fig1 = plt.figure(1, figsize=(8, 6*1.3))
fig1.set_facecolor('w')
sp1 = plt.subplot(3,1,1)
sp2 = plt.subplot(3,1,2, sharex=sp1)
sp3 = plt.subplot(3,1,3, sharex=sp1)


for ibuf, buf in enumerate(list_buffers_rec):
        ss = ch.buffer_2_beam(buf)
        #sp1.axvline(x=ss.slice_info['z_bin_center'], color='k', alpha=0.3, linestyle='--')
        sp1.axvspan(xmin=ss.slice_info['z_bin_left'], xmax=ss.slice_info['z_bin_right'],
            color={0:'green', 1:'orange'}[ibuf%2], alpha = 0.3)
        sp2.stem([ss.slice_info['z_bin_center']], [ss.slice_info['interact_with_EC']])
        sp3.stem([ss.slice_info['z_bin_center']], [ss.slice_info['i_slice']])
sp1.plot(thin_slice_set.z_centers, thin_slice_set.lambda_z(smoothen=False),'-', linewidth=2, color='k')

# check bunch merge
list_bunches_rec = []
for ibun, lbun in enumerate(list_bunch_buffers):
    list_sl = list(map(ch.buffer_2_beam, lbun))
    list_bunches_rec.append(st.merge_slices_into_bunch(list_sl))
    
fig2 = plt.figure(2, figsize=(8, 6*1.3))
fig2.set_facecolor('w')
spb1 = plt.subplot(3,1,1, sharex=sp1)

spb2 = plt.subplot(3,1,2, sharex=sp1)
spb3 = plt.subplot(3,1,3, sharex=sp1)
for ibun, bun in enumerate(list_bunches_rec):
        #spb1.axvline(x=bun.slice_info['z_bin_center'], color='k', alpha=0.3, linestyle='--')
        spb1.axvspan(xmin=bun.slice_info['z_bin_left'], xmax=bun.slice_info['z_bin_right'],
            color={0:'r', 1:'b'}[ibun%2], alpha = 0.3)
        spb2.stem([bun.slice_info['z_bin_center']], [bun.slice_info['interact_with_EC']])
        spb3.stem([bun.slice_info['z_bin_center']], [bun.slice_info['i_bunch']])
spb1.plot(thin_slice_set.z_centers, thin_slice_set.lambda_z(smoothen=False), linewidth=2, color='k')
        
for sp in [sp1, sp2, sp3, spb1, spb2, spb3]:
    sp.grid('on')

for sp in [sp1, spb1]:
    sp.set_ylabel('Beam profile [C/m]')
for sp in [sp2, spb2]:
    sp.set_ylabel('Interact with E.C.')

spb3.set_ylabel('Bunch number')
sp3.set_ylabel('Slice number')   

for sp in [sp3, spb3]:
    sp.set_xlabel('z [m]')

import pprint
# bunch to display
i_obs = 1 


bobs = list_bunches[-i_obs-1]
print('\n\nThe first bunch looks like:')
pprint.pprint(bobs.slice_info)
# Slice first bunch
if bobs.slice_info['slice_4_EC']:
    bobs_slices = st.slice_a_bunch(bobs, z_cut=z_cut, n_slices=n_slices)
    print('\n\nThe first slice of b0 looks like:')
    pprint.pprint(bobs_slices[-1].slice_info)
    print('\n\nThe last slice of b0 looks like:')
    pprint.pprint(bobs_slices[0].slice_info)
    print('\n\nThe mid slice of b0 looks like:')
    pprint.pprint(bobs_slices[n_slices//2].slice_info)

plt.show()



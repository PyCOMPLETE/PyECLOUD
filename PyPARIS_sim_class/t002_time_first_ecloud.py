# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

from . import Simulation as sim_mod

ring = sim_mod.get_serial_CPUring()

list_slice_objects = ring.pieces_to_be_treated # Head is the last element
list_machine_elements = ring.sim_content.mypart

# Truck the beam to the first ecloud (excluded)
for ee in list_machine_elements:
    if ee in ring.sim_content.my_list_eclouds:
        first_ecloud = ee
        break
    for ss in list_slice_objects[::-1]:
        ee.track(ss)

# Time first e-cloud
import time
import numpy as np
N_slices = len(list_slice_objects)

t_start = time.mktime(time.localtime())
for i_ss, ss in enumerate(list_slice_objects[::-1]):
    if np.mod(i_ss, 20)==0:
        print(("%d / %d"%(i_ss, N_slices)))
    first_ecloud.track(ss)
t_end = time.mktime(time.localtime())

print(('Track time %.2f s' % (t_end - t_start)))

import sys, os
sys.path.append('../../../')

import numpy as np
import matplotlib.pyplot as plt

import PyPARIS.myfilemanager as mfm

# first run 000, 001 or 002 script to generate footprint.h5 file

ob = mfm.object_with_arrays_and_scalar_from_h5('footprint.h5')

plt.plot(ob.qx_i,ob.qy_i,'ko')
plt.xlabel('Qx')
plt.ylabel('Qy')
plt.title('Scatter plot should resemble a triangle.')

print('Saving footprint.png.')
plt.savefig('footprint.png')

plt.show()

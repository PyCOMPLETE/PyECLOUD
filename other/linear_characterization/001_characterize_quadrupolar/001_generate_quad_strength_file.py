import os

import numpy as np
import matplotlib.pyplot as plt

import PyECLOUD.myfilemanager as mfm

plane = 'y'

# Load response data
ob = mfm.myloadmat_to_obj(f'./response_dc_{plane}.mat')

z_resp = ob.z_slices
r_resp_mat = ob.r_ideal
r_resp_mat[np.isnan(r_resp_mat)] = 0.
dpr_resp_mat = ob.dpr_slices_all_clouds
dpr_resp_mat[np.isnan(dpr_resp_mat)] = 0.

k_z = dpr_resp_mat/r_resp_mat

import scipy.io as sio
sio.savemat(f'linear_strength_{plane}.mat', {
    'z_slices': z_resp,
    'k_z_integrated': k_z})
plt.close('all')
fig100 = plt.figure(100)
ax101 = fig100.add_subplot(111)
ax101.plot(z_resp, k_z)

mask_fit = np.array(len(z_resp)*[True])
#if n_cut_fit > 0:
#    mask_fit[:n_cut_fit] = False
#    mask_fit[-n_cut_fit:] = False

p = np.polyfit(z_resp[mask_fit], k_z[mask_fit], deg=20)
ax101.plot(z_resp, np.polyval(p, z_resp))
## Temp for debug
#x_test = np.zeros(200)
#x_test[150] = 1
#dpx_test = np.dot(WW, x_test)
#figure(1000)
#plot(dpx_test)

plt.show()

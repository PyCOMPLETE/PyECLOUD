import pickle
import copy

import numpy as np
import matplotlib.pyplot as plt

import PyECLOUD.myfilemanager as mfm

n_sigma_ave = 2.5

ob = mfm.myloadmat_to_obj('field_map.mat')
ix_zero = np.argmin(np.abs(ob.xg))
iy_zero = np.argmin(np.abs(ob.yg))

mask_x = np.abs(ob.xg) <= n_sigma_ave*ob.sigma_x
mask_y = np.abs(ob.yg) <= n_sigma_ave*ob.sigma_y

Ex_lin = ob.Ex_L_map*0
Ey_lin = ob.Ey_L_map*0

k_x_list = []
k_y_list = []
for ss in range(len(ob.zg)):
    p_fit_x = np.polyfit(ob.xg[mask_x],
            ob.Ex_L_map[ss, mask_x, iy_zero], deg = 1)
    p_fit_y = np.polyfit(ob.yg[mask_y],
            ob.Ey_L_map[ss, ix_zero, mask_y], deg = 1)

    for iy in range(len(ob.yg)):
        Ex_lin[ss, :, iy] = p_fit_x[0]*ob.xg

    for ix in range(len(ob.xg)):
        Ey_lin[ss, ix :] = p_fit_y[0]*ob.yg

    k_x_list.append(p_fit_x[0])
    k_y_list.append(p_fit_y[0])

Ex_lin_ave = 0*Ex_lin
Ey_lin_ave = 0*Ex_lin
for ss in range(len(ob.zg)):
    for iy in range(len(ob.yg)):
        Ex_lin_ave[ss, :, iy] = np.mean(k_x_list)*ob.xg

    for ix in range(len(ob.xg)):
        Ey_lin_ave[ss, ix :] = np.mean(k_y_list)*ob.yg

Ex_nonlin = ob.Ex_L_map - Ex_lin
Ey_nonlin = ob.Ey_L_map - Ey_lin

Ex_nonlin_ave = 0*Ex_nonlin
Ey_nonlin_ave = 0*Ex_nonlin
for ss in range(len(ob.zg)):
    for iy in range(len(ob.yg)):
        Ex_nonlin_ave[ss, :, :] = np.mean(Ex_nonlin, axis=0)
    for ix in range(len(ob.xg)):
        Ey_nonlin_ave[ss, :, :] = np.mean(Ey_nonlin, axis=0)


dictlin = {kk: getattr(ob, kk) for kk in dir(ob) if not kk.startswith('__')}
dictlin['Ex_L_map'] = Ex_lin
dictlin['Ey_L_map'] = Ey_lin


dictlinave = {kk: getattr(ob, kk) for kk in dir(ob) if not kk.startswith('__')}
dictlinave['Ex_L_map'] = Ex_lin_ave
dictlinave['Ey_L_map'] = Ey_lin_ave

dictlinnoave = {kk: getattr(ob, kk) for kk in dir(ob) if not kk.startswith('__')}
dictlinnoave['Ex_L_map'] = Ex_lin - Ex_lin_ave
dictlinnoave['Ey_L_map'] = Ey_lin - Ey_lin_ave

dictnonlin = {kk: getattr(ob, kk) for kk in dir(ob) if not kk.startswith('__')}
dictnonlin['Ex_L_map'] = Ex_nonlin
dictnonlin['Ey_L_map'] = Ey_nonlin

dictnonlinave = {kk: getattr(ob, kk) for kk in dir(ob) if not kk.startswith('__')}
dictnonlinave['Ex_L_map'] = Ex_nonlin_ave
dictnonlinave['Ey_L_map'] = Ey_nonlin_ave

import scipy.io as sio
sio.savemat('field_map_lin.mat', dictlin)
sio.savemat('field_map_lin_ave.mat', dictlinave)
sio.savemat('field_map_lin_noave.mat', dictlinnoave)
sio.savemat('field_map_nonlin.mat', dictnonlin)
sio.savemat('field_map_nonlinave.mat', dictnonlinave)

i_slice = 100

plt.close('all')
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(211)
ax2 = fig1.add_subplot(212)

ax1.plot(ob.xg, ob.Ex_L_map[i_slice, :, iy_zero])
ax1.plot(ob.xg, Ex_lin[i_slice, :, iy_zero])
ax1.axvline(ob.sigma_x)
ax1.axvline(-ob.sigma_x)
ax1.set_ylim(2*np.array([-1, 1])*np.max(np.abs(Ex_lin[i_slice, mask_x, iy_zero])))


ax2.plot(ob.yg, ob.Ey_L_map[i_slice, ix_zero, :])
ax2.plot(ob.yg, Ey_lin[i_slice, ix_zero, :])
ax2.axvline(ob.sigma_y)
ax2.axvline(-ob.sigma_y)
ax2.set_ylim(2*np.array([-1, 1])*np.max(np.abs(Ey_lin[i_slice, ix_zero, mask_y])))

plt.show()

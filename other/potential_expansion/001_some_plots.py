import sys, os
sys.path.append(os.path.expanduser('../../../'))

import PyECLOUD.myfilemanager as mfm
import PyECLOUD.mystyle as ms

import numpy as np
import matplotlib.pyplot as plt

N_discard = 10

ob = mfm.myloadmat_to_obj('pinch_pic_data.mat')

ob.xg  =  ob.xg[N_discard:-N_discard]
ob.yg  =  ob.yg[N_discard:-N_discard]  
ob.rho =  ob.rho[:, N_discard:-N_discard, N_discard:-N_discard] 
ob.phi =  ob.phi[:, N_discard:-N_discard, N_discard:-N_discard] 
ob.Ex  =  ob.Ex[:, N_discard:-N_discard, N_discard:-N_discard] 
ob.Ey  =  ob.Ey[:, N_discard:-N_discard, N_discard:-N_discard] 

x_obs = 0.
ix_obs = np.argmin(np.abs(ob.xg - x_obs))
y_obs = 0.
iy_obs = np.argmin(np.abs(ob.yg - y_obs))
z_obs = 0.2
iz_obs = np.argmin(np.abs(ob.zg - z_obs))


plt.close('all')
ms.mystyle_arial()
fig1 = plt.figure()
fig1.set_facecolor('w')
ax1 = fig1.add_subplot(3,1,1)
ax2 = fig1.add_subplot(3,1,2, sharex=ax1)
ax3 = fig1.add_subplot(3,1,3, sharex=ax1)

ax1.pcolormesh(ob.zg, ob.yg, ob.rho[:,ix_obs,:].T)
ax2.pcolormesh(ob.zg, ob.yg, ob.Ey[:,ix_obs,:].T)
ax3.pcolormesh(ob.zg, ob.yg, ob.phi[:,ix_obs,:].T)

fig2 = plt.figure(2, figsize=(8*1.5, 6*1.5))
fig2.set_facecolor('w')
axc1 = plt.subplot2grid(shape=(3, 3), loc=(0,0), colspan=2, fig=fig2)
axc2 = plt.subplot2grid(shape=(3, 3), loc=(1,0), colspan=2, fig=fig2, sharex=axc1)
axc3 = plt.subplot2grid(shape=(3, 3), loc=(2,0), colspan=2, fig=fig2, sharex=axc1)

axc1.plot(ob.yg, ob.rho[iz_obs, ix_obs, :].T, '.-')
axc1.plot(ob.xg, ob.rho[iz_obs, :, iy_obs].T, '.-r')
axc2.plot(ob.yg, ob.Ey[iz_obs, ix_obs,:].T, '.-')
axc2.plot(ob.xg, ob.Ex[iz_obs, :, iy_obs].T, '.-r')
axc3.plot(ob.yg, ob.phi[iz_obs, ix_obs,:].T, '.-')
axc3.plot(ob.xg, ob.phi[iz_obs, :, iy_obs].T, '.-r')

axc1.set_ylim(top=0.)

for ax in [axc1, axc2, axc3]:
    ax.ticklabel_format(style='sci', scilimits=(0,0),axis='y')
    ax.ticklabel_format(style='sci', scilimits=(0,0),axis='x')

fig2.subplots_adjust(hspace=.4)

plt.show()


import sys, os
sys.path.append(os.path.expanduser('../../../'))

import PyECLOUD.myfilemanager as mfm
import PyECLOUD.mystyle as ms

import numpy as np
import matplotlib.pyplot as plt

from scipy.constants import e as qe

N_discard = 10

ob = mfm.myloadmat_to_obj('pinch_pic_data.mat')

ob.xg  =  ob.xg[N_discard:-N_discard]
ob.yg  =  ob.yg[N_discard:-N_discard]  
ob.rho =  ob.rho[:, N_discard:-N_discard, N_discard:-N_discard] 
ob.phi =  ob.phi[:, N_discard:-N_discard, N_discard:-N_discard] 
ob.Ex  =  ob.Ex[:, N_discard:-N_discard, N_discard:-N_discard] 
ob.Ey  =  ob.Ey[:, N_discard:-N_discard, N_discard:-N_discard] 

ob.rho /= -qe

x_obs = 0.
ix_obs = np.argmin(np.abs(ob.xg - x_obs))
y_obs = 0.
iy_obs = np.argmin(np.abs(ob.yg - y_obs))

plt.close('all')
ms.mystyle_arial()
fig1 = plt.figure(1, figsize=(8, 6*1.5))
fig1.set_facecolor('w')
ax1 = fig1.add_subplot(3,1,1)
ax2 = fig1.add_subplot(3,1,2, sharex=ax1)
ax3 = fig1.add_subplot(3,1,3, sharex=ax1)

mbl = ax1.pcolormesh(ob.zg, ob.yg, ob.rho[:,ix_obs,:].T)
plt.colorbar(mappable=mbl, ax=ax1, aspect=5)
mbl = ax2.pcolormesh(ob.zg, ob.yg, ob.Ey[:,ix_obs,:].T)
plt.colorbar(mappable=mbl, ax=ax2, aspect=5, format='%.1e')
mbl = ax3.pcolormesh(ob.zg, ob.yg, ob.phi[:,ix_obs,:].T)
plt.colorbar(mappable=mbl, ax=ax3, aspect=5)

fig1.subplots_adjust(hspace=.4)



fig2 = plt.figure(2, figsize=(8*1.5, 6*1.5))

for i_frame, z_obs in enumerate(ob.zg[::-1]):

    fig2.clear()
    fig2.set_facecolor('w')
    
    iz_obs = np.argmin(np.abs(ob.zg - z_obs))
    
    axc1 = plt.subplot2grid(shape=(3, 3), loc=(0,0), colspan=2, fig=fig2)
    axc2 = plt.subplot2grid(shape=(3, 3), loc=(1,0), colspan=2, fig=fig2, sharex=axc1)
    axc3 = plt.subplot2grid(shape=(3, 3), loc=(2,0), colspan=2, fig=fig2, sharex=axc1)
    
    axc1.plot(ob.yg, ob.rho[iz_obs, ix_obs, :].T, '.-')
    axc1.plot(ob.xg, ob.rho[iz_obs, :, iy_obs].T, '.-r')
    axc1.set_ylim(np.min(ob.rho), np.max(ob.rho))
    
    axc2.plot(ob.yg, ob.Ey[iz_obs, ix_obs,:].T, '.-')
    axc2.plot(ob.xg, ob.Ex[iz_obs, :, iy_obs].T, '.-r')
    axc2.set_ylim(np.min(ob.Ex), np.max(ob.Ex))
    
    axc3.plot(ob.yg, ob.phi[iz_obs, ix_obs,:].T, '.-')
    axc3.plot(ob.xg, ob.phi[iz_obs, :, iy_obs].T, '.-r')
    #axc3.set_ylim(np.min(ob.phi), np.max(ob.phi))
    
    for ax in [axc1, axc2, axc3]:
        ax.axvline(x=ob.sigma_y_beam, color='b', linestyle='--')
        ax.axvline(x=-ob.sigma_y_beam, color='b', linestyle='--')
    
    axd1 = plt.subplot2grid(shape=(3, 3), loc=(0,2), colspan=1, fig=fig2)
    axd1.pcolormesh(ob.xg, ob.yg, ob.rho[iz_obs, :, :].T,
            vmin=np.min(ob.rho), vmax=np.max(ob.rho))
    axd1.axis('equal')
    
    
    for ax in [axc1, axc2, axc3, ax1, ax2, ax3, axd1]:
        ax.ticklabel_format(style='sci', scilimits=(0,0),axis='y')
        ax.ticklabel_format(style='sci', scilimits=(0,0),axis='x')
    
    fig2.subplots_adjust(hspace=.4, wspace=.32)
    
    fig2.suptitle('z = %.2e (%.2e sigmaz)'%(ob.zg[iz_obs], 
        ob.zg[iz_obs]/ob.sigma_z_beam))
    
    fig2.savefig('temp/frame_%03d.png'%i_frame, dpi=150)

plt.show()


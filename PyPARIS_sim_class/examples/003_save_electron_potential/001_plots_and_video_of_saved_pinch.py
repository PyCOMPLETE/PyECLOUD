import sys, os
sys.path.append('../../../')

import shutil
import PyECLOUD.myfilemanager as mfm
import PyECLOUD.mystyle as ms

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from scipy.constants import e as qe

N_discard = 10
tag = 'movie_ecloud'
gen_movie = True 
z_single_frame_interactive = 1e-2 

# ======================

folder_movie = 'temp_'+tag+'_frames'
os.mkdir(folder_movie)

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

z_movie = ob.zg[::-1]


fig2 = plt.figure(2, figsize=(8*1.5, 6*1.5))

for i_frame, z_obs in enumerate(list(z_movie) + [z_single_frame_interactive]):
    
    print(('Frame %d/%d'%(i_frame, len(z_movie))))

    if not gen_movie and i_frame < len(z_movie):
        continue

    fig2.clear()
    fig2.set_facecolor('w')
    
    iz_obs = np.argmin(np.abs(ob.zg - z_obs))
    
    axc1 = plt.subplot2grid(shape=(3, 3), loc=(0,0), colspan=2, fig=fig2)
    axc2 = plt.subplot2grid(shape=(3, 3), loc=(1,0), colspan=2, fig=fig2, sharex=axc1)
    axc3 = plt.subplot2grid(shape=(3, 3), loc=(2,0), colspan=2, fig=fig2, sharex=axc1)
    
    axc1.plot(ob.yg, ob.rho[iz_obs, ix_obs, :].T, '.-', label='x=0')
    axc1.plot(ob.xg, ob.rho[iz_obs, :, iy_obs].T, '.-r', label='y=0')
    axc1.set_ylim(np.min(ob.rho), np.max(ob.rho))
    axc1.set_ylabel('rho [e-/m^3]')
    axc1.legend(loc='upper right', prop={'size':14})
    
    axc2.plot(ob.yg, ob.Ey[iz_obs, ix_obs,:].T, '.-')
    axc2.plot(ob.xg, ob.Ex[iz_obs, :, iy_obs].T, '.-r')
    axc2.set_ylim(np.min(ob.Ex), np.max(ob.Ex))
    axc2.set_ylabel('Electric field [V/m]')

    axc3.plot(ob.yg, ob.phi[iz_obs, ix_obs,:].T, '.-')
    axc3.plot(ob.xg, ob.phi[iz_obs, :, iy_obs].T, '.-r')
    #axc3.set_ylim(np.min(ob.phi), np.max(ob.phi))
    axc3.set_ylabel('phi [V]')
    axc3.set_xlabel('Position [m]')
    
    for ax in [axc1, axc2, axc3]:
        ax.axvline(x=ob.sigma_y_beam, color='b', linestyle='--')
        ax.axvline(x=-ob.sigma_y_beam, color='b', linestyle='--')
    
    axd1 = plt.subplot2grid(shape=(3, 3), loc=(0,2), colspan=1, fig=fig2)
    mbl = axd1.pcolormesh(ob.xg*1e3, ob.yg*1e3, 
            np.log10(ob.rho[iz_obs, :, :].T),
            vmin=np.log10(np.max(ob.rho))-3., 
            vmax=np.log10(np.max(ob.rho)))
    axd1.axis('equal')
    plt.colorbar(mappable=mbl, ax=axd1, aspect=20)
    axd1.yaxis.set_major_locator(MaxNLocator(5))
    axd1.xaxis.set_major_locator(MaxNLocator(5))
    axd1.set_title('log10(rho)')
    
    axd3 = plt.subplot2grid(shape=(3, 3), loc=(2,2), colspan=1, fig=fig2,
            sharex=axd1)
    mbl = axd3.pcolormesh(ob.xg*1e3, ob.yg*1e3, 
            ob.phi[iz_obs, :, :].T)
    axd3.axis('equal')
    plt.colorbar(mappable=mbl, ax=axd3, aspect=20)
    # axd3.yaxis.set_major_locator(MaxNLocator(5))
    # axd3.xaxis.set_major_locator(MaxNLocator(5))
    
    for ax in [axc1, axc2, axc3, ax1, ax2, ax3, axd1]:
        ax.ticklabel_format(style='sci', scilimits=(0,0),axis='y')
        ax.ticklabel_format(style='sci', scilimits=(0,0),axis='x')
    
    fig2.subplots_adjust(hspace=.4, wspace=.32, right=.96)
    
    fig2.suptitle('z = %.2e (%.2f sigmaz)'%(ob.zg[iz_obs], 
        ob.zg[iz_obs]/ob.sigma_z_beam))
    
    if gen_movie and i_frame < len(z_movie):
        fig2.savefig(folder_movie+'/frame_%03d.png'%i_frame, dpi=150)


if gen_movie:
    os.system(' '.join([
        'ffmpeg',
        '-i %s'%folder_movie+'/frame_%03d.png',
        '-c:v libx264 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2,setpts=4.*PTS"',
        '-profile:v high -level:v 4.0 -pix_fmt yuv420p -crf 22',
        '-codec:a aac movie_%s.mp4'%tag])) 
   
    shutil.rmtree(folder_movie)

plt.show()


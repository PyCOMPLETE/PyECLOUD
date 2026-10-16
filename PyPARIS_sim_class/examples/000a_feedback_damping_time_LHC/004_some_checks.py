import sys, os
sys.path.append('../../../')


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

import PyPARIS.myfilemanager as mfm

import PyECLOUD.mystyle as ms

tag = 'noecloud'

flag_check_damp_time = True
tau_damp_x = 50.
tau_damp_y = 100.

flag_check_Qs = True
Q_s = np.nan

#ob = mfm.monitorh5_to_obj('bunch_evolution_00.h5')
ob = mfm.monitorh5list_to_obj(
        ['bunch_evolution_%02d.h5'%ii for ii in range(3)])


plt.close('all')

# Plot transverse positions
fig1 = plt.figure(1, figsize=(8,6*1.2))
fig1.set_facecolor('w')
ms.mystyle_arial(fontsz=16)

axx = fig1.add_subplot(2,1,1)
axy = fig1.add_subplot(2,1,2, sharex=axx)

axx.plot(ob.mean_x)
axy.plot(ob.mean_y)

if flag_check_damp_time:
    turn_num = np.arange(0, len(ob.mean_x), dtype=np.float)
    ix_max = np.argmax(ob.mean_x)
    iy_max = np.argmax(ob.mean_y)

    axx.plot(ob.mean_x[ix_max]*np.exp(-(turn_num-ix_max)/tau_damp_x),
            linewidth=2, color='red', linestyle='--',
            label=r'Damping time = %.0f turns'%tau_damp_x)
    axy.plot(ob.mean_y[iy_max]*np.exp(-(turn_num-iy_max)/tau_damp_y),
            linewidth=2, color='red', linestyle='--',
            label=r'Damping time = %.0f turns'%tau_damp_y)

axx.legend(prop={'size':14})
axy.legend(prop={'size':14})

axx.set_ylabel('x [m]')
axy.set_ylabel('y [m]')
axy.set_xlabel('Turn')

# Plot transverse spectra
fig2 = plt.figure(2, figsize=(8,6*1.2))
fig2.set_facecolor('w')

axfx = fig2.add_subplot(2,1,1)
axfy = fig2.add_subplot(2,1,2, sharex=axfx)

spectx = np.abs(np.fft.rfft(ob.mean_x))
specty = np.abs(np.fft.rfft(ob.mean_y))
freq = np.fft.rfftfreq(len(ob.mean_x))

axfx.plot(freq, spectx)
axfy.plot(freq, specty)

# Check longitudinal plane
fig3 = plt.figure(3, figsize=(8,6*1.2))
fig3.set_facecolor('w')
axz = fig3.add_subplot(2,1,1, sharex=axx)
axfz = fig3.add_subplot(2,1,2)

axz.plot(ob.mean_z[:-10])
spectz = np.abs(np.fft.rfft(ob.mean_z[:-10]))
spectz[0] = 0. # I am non interested in the average
freqz = np.fft.rfftfreq(len(ob.mean_x[:-10]))
axfz.plot(freqz, spectz)
axfz.axvline(x=Q_s)

for ax in [axx, axy, axfx, axfy, axz, axfz]:
    ax.ticklabel_format(style='sci', scilimits=(0,0),axis='y')

for ax in [axx, axy, axfx, axfy, axz, axfz]:
    ax.yaxis.set_major_locator(MaxNLocator(5))

for ax in [axx, axy, axfx, axfy, axz, axfz]:
    ax.grid(True)

for fig in [fig1, fig2, fig3]:
    fig.suptitle('Bunch')
    fig.subplots_adjust(
            top=0.885,
            bottom=0.1,
            left=0.125,
            right=0.9,
            hspace=0.345,
            wspace=0.2)

plt.show()

import sys, os
sys.path.append('../../../')


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

import PyPARIS.myfilemanager as mfm

import PyECLOUD.mystyle as ms

tag = 'noecloud'
i_bunch = 0

flag_check_damp_time = True
tau_damp_x = 20.
tau_damp_y = 200000.

flag_check_Qs = True
Q_s = 5.664e-03

ob = mfm.myloadmat_to_obj(tag+'_matrices.mat')

plt.close('all')

# Plot transverse positions
fig1 = plt.figure(1, figsize=(8,6*1.2))
fig1.set_facecolor('w')
ms.mystyle_arial(fontsz=16)

axx = fig1.add_subplot(2,1,1)
axy = fig1.add_subplot(2,1,2, sharex=axx)

axx.plot(ob.mean_x[:, i_bunch])
axy.plot(ob.mean_y[:, i_bunch])

if flag_check_damp_time:
    turn_num = np.arange(0, len(ob.mean_x[:, i_bunch]), dtype=np.float)
    axx.plot(ob.mean_x[0, i_bunch]*np.exp(-turn_num/tau_damp_x),
            linewidth=2, color='red', linestyle='--',
            label=r'Damping time = %.0f turns'%tau_damp_x)
    axy.plot(ob.mean_y[0, i_bunch]*np.exp(-turn_num/tau_damp_y),
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

spectx = np.abs(np.fft.rfft(ob.mean_x[:, i_bunch]))
specty = np.abs(np.fft.rfft(ob.mean_y[:, i_bunch]))
freq = np.fft.rfftfreq(len(ob.mean_x[:, i_bunch]))

axfx.plot(freq, spectx)
axfy.plot(freq, specty)

# Check longitudinal plane
fig3 = plt.figure(3, figsize=(8,6*1.2))
fig3.set_facecolor('w')
axz = fig3.add_subplot(2,1,1, sharex=axx)
axfz = fig3.add_subplot(2,1,2)

axz.plot(ob.mean_z[:-10, i_bunch])
spectz = np.abs(np.fft.rfft(ob.mean_z[:-10, i_bunch]))
spectz[0] = 0. # I am non interested in the average
freqz = np.fft.rfftfreq(len(ob.mean_x[:-10, i_bunch]))
axfz.plot(freqz, spectz)
axfz.axvline(x=Q_s)

for ax in [axx, axy, axfx, axfy, axz, axfz]:
    ax.ticklabel_format(style='sci', scilimits=(0,0),axis='y')

for ax in [axx, axy, axfx, axfy, axz, axfz]:
    ax.yaxis.set_major_locator(MaxNLocator(5))

for ax in [axx, axy, axfx, axfy, axz, axfz]:
    ax.grid(True)

for fig in [fig1, fig2, fig3]:
    fig.suptitle('Bunch %d'%i_bunch)
    fig.subplots_adjust(
            top=0.885,
            bottom=0.1,
            left=0.125,
            right=0.9,
            hspace=0.345,
            wspace=0.2)

plt.show()

import sys, os
BIN = os.path.expanduser("../../../")
sys.path.append(BIN)

import argparse
import pylab as pl
import numpy as np
#from colorsys import hsv_to_rgb
import os
import PyECLOUD.myloadmat_to_obj as mlm
import matplotlib.gridspec as gridspec
import PyECLOUD.mystyle as ms

pl.close('all')

sim_folder = './'


dict_ref = mlm.myloadmat('../tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_multigrid/' + '/Pyecltest_angle3D_ref.mat') # load dictionary of the reference simulation
dict_curr0 = mlm.myloadmat('./test_saving__iter0.mat')   # load dictionary of the current simulation
dict_curr1 = mlm.myloadmat('./test_saving__iter1.mat')   # load dictionary of the current simulation

pl.close('all')
myfontsz = 16
ms.mystyle_arial(fontsz=myfontsz)

fig1 = pl.figure(1)
fig1.set_facecolor('w')
ax1 = fig1.add_subplot(1,1,1)
ax1.plot(dict_ref['t'], dict_ref['Nel_timep'], '-', linewidth=2,
        label="Buildup simulation")
ax1.plot(dict_curr0['t'] - 10e-9 - dict_curr0['t'][0], dict_curr0['Nel_timep'],
        'r-', linewidth=2, label='PyHT beam, turn 1')
ax1.plot(dict_curr1['t'] - 10e-9 - dict_curr1['t'][0], dict_curr1['Nel_timep'], 
        'g-', linewidth=2, label='PyHT beam, turn 2')
ax1.ticklabel_format(style='sci', scilimits=(0, 0), axis='x')
ax1.grid(True)
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Number of electrons [m^-1]')
ax1.legend(prop={'size':14}, loc='upper left')

fig1.subplots_adjust(bottom=.12)


fig2 = pl.figure(2)
fig2.set_facecolor('w')
ax2 = fig2.add_subplot(1,1,1)

ax2.plot(np.sum(dict_ref['En_hist'], axis=1), '.-', linewidth=2,
        markersize=10, label="Buildup simulation")
ax2.plot(np.sum(dict_curr0['En_hist'], axis=1), '.r-', linewidth=2,
        markersize=10, label='PyHT beam, turn 1')
ax2.plot(np.sum(dict_curr1['En_hist'], axis=1), '.g-', linewidth=2,
        markersize=10, label='PyHT beam, turn 2')

ax2.set_xlabel('Bunch passage')
ax2.set_ylabel('Heat load eV/bunch')
ax2.legend(prop={'size':14}, loc='upper left')
ax2.grid(True)
dict_curr = dict_curr1

#~ pl.figure(1001)
#~ pl.semilogy(dict_ref['t'], dict_ref['En_kin_eV_time'], '.-')
#~ pl.semilogy(dict_curr['t']-10e-9, dict_curr['En_kin_eV_time'], '.r-')

#~ pl.figure(1002)
#~ pl.plot(dict_ref['t'][:-1], np.diff(dict_ref['t'])/1e-12, '.-')
#~ pl.plot(dict_curr['t'][:-1]-10e-9, np.diff(dict_curr['t'])/1e-12, '.r-')

pl.figure(1003)
pl.plot(dict_ref['N_mp_pass'], '.-')
pl.plot(dict_curr['N_mp_pass'], '.r-')

pl.figure(1004)
i_pass = 5
pl.plot(dict_ref['nel_hist'][i_pass, :], '.-')
pl.plot(dict_curr['nel_hist'][i_pass, :], '.r-')

#~ pl.figure(1005)
#~ pl.plot(dict_ref['En_hist'][i_pass, :], '.-')
#~ pl.plot(dict_curr['En_hist'][i_pass, :], '.r-')

#~ pl.figure(1006)
#~ pl.plot(dict_ref['t'][:-1], dict_ref['En_imp_eV_time'][:-1]/np.diff(dict_ref['t']), '.-')
#~ pl.plot(dict_curr['t'][:-1]-10e-9, dict_curr['En_imp_eV_time'][:-1]/np.diff(dict_curr['t']), '.r-')


#~ pl.figure(1007)
#~ pl.plot(dict_ref['energ_eV_impact_hist'][i_pass, :], '.-')
#~ pl.plot(dict_curr['energ_eV_impact_hist'][i_pass, :], '.r-')

#~ pl.figure(1008)
#~ pl.semilogy(dict_ref['t'], dict_ref['lam_t_array'], '.-')
#~ pl.semilogy(dict_curr['t']-10e-9, dict_curr['lam_t_array'], '.r-')

pl.figure(1009)
pl.plot(np.sum(dict_ref['nel_hist'], axis=1), '.-')
pl.plot(np.sum(dict_curr['nel_hist'], axis=1), '.r-')


pl.figure(1010)
pl.plot(dict_ref['N_mp_ref_pass'], '.-')
pl.plot(dict_curr['N_mp_ref_pass'], '.r-')


pl.show()
pl.show()




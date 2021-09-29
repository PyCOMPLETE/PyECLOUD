import sys
sys.path.append('../../../')


import pylab as pl
import numpy as np
import seaborn as sns
import time

import PyECLOUD.myfilemanager as mfm

pyhdtl = mfm.monitorh5_to_dict('bunch_evolution_A44_156b_26ips_10turns_5.00nTorr.h5',
        key='Slices')

N_turns = 10

pl.close('all')
pl.ion()

sns.set_context('talk', font_scale=1.4, rc={'lines.linewidth': 1.5})
sns.set_style('whitegrid', {'grid.linestyle': ':', 'axes.edgecolor': '0.5', 'axes.linewidth': 1.2, 'legend.frameon': False})
mksize = 7
fsize = 12


pl.figure(100, figsize=(8, 10))
sp3 = pl.subplot(2, 1, 1)
pl.plot(pyhdtl['mean_x'][::-1, 0], '.-', markersize=mksize)
pl.gca().ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
pl.ylabel('Horizontal centroid position\n [m]')
pl.xlabel('Bunch nr')
pl.xlim(0, 156)
pl.title(('After %d turns')%(0))

sp4 = pl.subplot(2, 1, 2, sharex=sp3, sharey=sp3)
pl.plot(pyhdtl['mean_y'][::-1, 0], '.-', markersize=mksize)
pl.gca().ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
pl.ylabel('Vertical centroid position\n [m]')
pl.xlabel('Bunch nr')
pl.xlim(0, 156)
pl.ylim(-3e-6, 3e-6)
pl.subplots_adjust(left=0.21, hspace=0.3)


for ii in range(N_turns + 1):

	pl.clf()
	sp3 = pl.subplot(2, 1, 1)
	pl.plot(pyhdtl['mean_x'][::-1, ii], '.-', markersize=mksize)
	pl.gca().ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
	pl.ylabel('Horizontal centroid position\n [m]')
	pl.xlabel('Bunch nr')
	pl.xlim(0, 156)
	pl.title(('After %d turns')%(ii))

	sp4 = pl.subplot(2, 1, 2, sharex=sp3, sharey=sp3)
	pl.plot(pyhdtl['mean_y'][::-1, ii], '.-', markersize=mksize)
	pl.gca().ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
	pl.ylabel('Vertical centroid position\n [m]')
	pl.xlabel('Bunch nr')
	pl.xlim(0, 156)
	pl.ylim(-3e-6, 3e-6)
	pl.subplots_adjust(left=0.21, hspace=0.3)

	pl.pause(1.5)



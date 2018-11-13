import pylab as pl
import numpy as np
import seaborn as sns
import myfilemanager as mfm
import pickle


pyhdtl = mfm.bunchh5_to_dict('bunch_evolution_A44_156b_260ips_1turns_20.00nTorr.h5')
lattice_file = 'CLIC_DR_n260_optics.pkl'


with open(lattice_file, 'r') as fid:
	optics = pickle.load(fid)

pl.close('all')

sns.set_context('talk', font_scale=1.2, rc={'lines.linewidth': 1.5})
sns.set_style('whitegrid', {'grid.linestyle': ':', 'axes.edgecolor': '0.5', 'axes.linewidth': 1.2})
mksize = 7


ip = 258
pl.figure(figsize=(7.5, 10))
sp3 = pl.subplot(2, 1, 1)
pl.plot(pyhdtl['mean_x'][::-1, ip], '.-', markersize=mksize)
pl.gca().ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
pl.ylabel('Horizontal centroid position\n [m]')
pl.xlabel('Bunch nr')
pl.xlim(0, 156)
pl.title(('Bunch train at IP %d')%(ip))

sp4 = pl.subplot(2, 1, 2, sharex=sp3, sharey=sp3)
pl.plot(pyhdtl['mean_y'][::-1, ip], '.-', markersize=mksize)
pl.gca().ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
pl.ylabel('Vertical centroid position\n [m]')
pl.xlabel('Bunch nr')
pl.xlim(0, 156)
pl.subplots_adjust(left=0.21, hspace=0.3)


i_bunch = 156
n_bunch = len(pyhdtl['mean_x'][:, 0])
pl.figure(figsize=(7.5, 10))
sp1 = pl.subplot(2, 1, 1)
pl.plot(optics['s'][1:-1], pyhdtl['mean_x'][n_bunch - i_bunch, :-1] / pyhdtl['sigma_x'][n_bunch - i_bunch, :-1], '.-', markersize=mksize, linewidth=1.6)
pl.ylabel('Bunch centroid position\n $<x>/\sigma_x$')
pl.xlabel('s [m]')
pl.title(('Bunch %d')%(i_bunch))

sp2 = pl.subplot(2, 1, 2, sharex=sp1, sharey=sp1)
pl.plot(optics['s'][1:-1], pyhdtl['mean_y'][n_bunch - i_bunch, :-1] / pyhdtl['sigma_y'][n_bunch - i_bunch, :-1], '.-', markersize=mksize, linewidth=1.6)
pl.ylabel('Bunch centroid position\n $<y>/\sigma_y$')
pl.xlabel('s [m]')
pl.xlim(0, 428)
pl.subplots_adjust(left=0.21, hspace=0.3)


pl.show()

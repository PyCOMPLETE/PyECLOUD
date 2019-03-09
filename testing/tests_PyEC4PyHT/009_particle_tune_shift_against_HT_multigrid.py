import sys, os
sys.path.append(os.path.expanduser('../../../'))
sys.path.append(os.path.expanduser('../../../PyHEADTAIL/'))


import numpy as np
import pylab as pl
import PyECLOUD.mystyle as ms

n_segments = 5
B_multip = [0.5]
n_record = 1000
n_part_per_turn = 5000


# define machine for PyHEADTAIL
from PyHEADTAIL.particles.slicing import UniformBinSlicer
from machines_for_testing  import shortSPS
machine = shortSPS(n_segments=n_segments, machine_configuration='Q20-injection-like')

# remove synchrotron motion
machine.one_turn_map.remove(machine.longitudinal_map)


def dict_of_arrays_and_scalar_from_h5(filename):
	import h5py
	with h5py.File(filename, 'r') as fid:
		f_dict = {}
		for kk in fid.keys():
			f_dict[kk] = np.array(fid[kk]).copy()
			if f_dict[kk].shape == ():
				f_dict[kk] = f_dict[kk].tolist()
	return f_dict

dict_HT = dict_of_arrays_and_scalar_from_h5('footprint_HT.h5ref')

qx_ht = dict_HT['qx_ht']
qy_ht = dict_HT['qy_ht']
qx_centroid_ht = dict_HT['qx_centroid_ht']
qy_centroid_ht = dict_HT['qy_centroid_ht']

# compute sigma x and y
epsn_x = 2.5e-6
epsn_y = 2.5e-6

inj_optics = machine.transverse_map.get_injection_optics()
sigma_x = np.sqrt(inj_optics['beta_x'] * epsn_x / machine.betagamma)
sigma_y = np.sqrt(inj_optics['beta_y'] * epsn_y / machine.betagamma)


# generate a bunch
bunch = machine.generate_6D_Gaussian_bunch(n_macroparticles=300000, intensity=1.15e11, epsn_x=epsn_x, epsn_y=epsn_y, sigma_z=0.2)

# replace first particles with HEADTAIL ones
bunch.x[:n_part_per_turn] = dict_HT['x0_HT']
bunch.xp[:n_part_per_turn] = dict_HT['xp0_HT']
bunch.y[:n_part_per_turn] = dict_HT['y0_HT']
bunch.yp[:n_part_per_turn] = dict_HT['yp0_HT']
bunch.z[:n_part_per_turn] = dict_HT['z0_HT']
bunch.dp[:n_part_per_turn] = dict_HT['dp0_HT']
n_turns = dict_HT['n_turns']


# define apertures and Dh_sc to simulate headtail conditions
x_aper = 20 * sigma_x
y_aper = 20 * sigma_y
Dh_sc = 2 * x_aper / 128 / 2.

# ecloud
import PyECLOUD.PyEC4PyHT as PyEC4PyHT
from PyHEADTAIL.particles.slicing import UniformBinSlicer
slicer = UniformBinSlicer(n_slices=64, z_cuts=(-3 * bunch.sigma_z(), 3 * bunch.sigma_z()))


init_unif_edens_flag = 1
init_unif_edens = 1e11
N_MP_ele_init = 100000
N_mp_max = N_MP_ele_init * 4.

nel_mp_ref_0 = init_unif_edens * 4 * x_aper * y_aper / N_MP_ele_init
Dh_ext = 2e-3

ecloud = PyEC4PyHT.Ecloud(L_ecloud=machine.circumference / machine.transverse_map.n_segments, slicer=slicer,
                          Dt_ref=25e-12, pyecl_input_folder='./drift_sim',
                          x_aper=x_aper, y_aper=y_aper,
                          init_unif_edens_flag=init_unif_edens_flag,
                          init_unif_edens=init_unif_edens,
                          N_mp_max=N_mp_max,
                          nel_mp_ref_0=nel_mp_ref_0,
                          B_multip=B_multip,
                          PyPICmode='ShortleyWeller_WithTelescopicGrids',
                          f_telescope=0.3,
                          Dh_sc=Dh_ext,
                          target_grid={'x_min_target': -4 * bunch.sigma_x(), 'x_max_target': 4 * bunch.sigma_x(),
                                       'y_min_target': -4 * bunch.sigma_y(), 'y_max_target': 4 * bunch.sigma_y(),
                                       'Dh_target': Dh_sc},
                          N_nodes_discard=10.,
                          N_min_Dh_main=10,)

# generate a bunch
bunch_for_map = machine.generate_6D_Gaussian_bunch(n_macroparticles=500000,
                                                   intensity=1.15e11, epsn_x=epsn_x, epsn_y=epsn_y, sigma_z=0.2)
ecloud.track_once_and_replace_with_recorded_field_map(bunch_for_map)

machine.install_after_each_transverse_segment(ecloud)


# prepare storage for particles cohordinates
x_i = np.empty((n_record, n_turns))
xp_i = np.empty((n_record, n_turns))
y_i = np.empty((n_record, n_turns))
yp_i = np.empty((n_record, n_turns))

# track and store
for i in range(n_turns):
    machine.track(bunch)  # , verbose=True)

    print 'Turn', i
    sys.stdout.flush()

    x_i[:, i] = bunch.x[:n_record]
    xp_i[:, i] = bunch.xp[:n_record]
    y_i[:, i] = bunch.y[:n_record]
    yp_i[:, i] = bunch.yp[:n_record]
print '\nDONE'

from tune_analysis import tune_analysis
qx_i, qy_i, qx_centroid, qy_centroid = tune_analysis(x_i, xp_i, y_i, yp_i)

pl.close('all')
ms.mystyle_arial(fontsz=14)
pl.figure(1)
sp1 = pl.subplot(2, 1, 1)
pl.plot(np.mean(x_i, axis=0), '.-b', markersize=5, linewidth=2, label='PyHT')
pl.ylabel('<x>')
pl.grid('on')
ms.sciy()
pl.legend(prop={'size': 14})
pl.subplot(2, 1, 2, sharex=sp1)
pl.plot(np.mean(y_i, axis=0), '.-b', markersize=5, linewidth=2, label='PyHT')
pl.xlabel('Turn'); pl.ylabel('<y>')
pl.grid('on')
ms.sciy()
#pl.savefig(filename.split('_prb.dat')[0]+'_centroids.png', dpi=200)


pl.figure(2)
pl.plot(np.abs(qx_i), np.abs(qy_i), '.', label='PyHT', markersize=3)
pl.plot(np.abs(qx_ht), np.abs(qy_ht), '.r', label='HT', markersize=3)
pl.plot([np.modf(machine.Q_x)[0]], [np.modf(machine.Q_y)[0]], 'go')
pl.xlabel('$Q_x$'); pl.ylabel('$Q_y$')
pl.legend(prop={'size': 14})
pl.grid('on')
pl.axis('equal')
#pl.savefig(filename.split('_prb.dat')[0]+'_footprint.png', dpi=200)


pl.figure(3)
pl.subplot(2, 1, 1)
pl.plot(bunch.z[:n_record], np.abs(qx_i) - np.modf(machine.Q_x)[0], '.', markersize=3, label='PyHT')
pl.plot(dict_HT['z0_HT'][:n_record], np.abs(qx_ht) - np.modf(machine.Q_x)[0], '.r', markersize=3, label='HT')
pl.ylabel('$\Delta Q_x$')
pl.grid('on')
pl.legend(prop={'size': 14})
pl.subplot(2, 1, 2)
pl.plot(bunch.z[:n_record], np.abs(qy_i) - np.modf(machine.Q_x)[0], '.', markersize=3)
pl.plot(dict_HT['z0_HT'][:n_record], np.abs(qy_ht) - np.modf(machine.Q_x)[0], '.r', markersize=3)
pl.ylabel('$\Delta Q_y$')
pl.xlabel('z [m]')
pl.grid('on')
#pl.savefig(filename.split('_prb.dat')[0]+'_Q_vs_z.png', dpi=200)
pl.show()

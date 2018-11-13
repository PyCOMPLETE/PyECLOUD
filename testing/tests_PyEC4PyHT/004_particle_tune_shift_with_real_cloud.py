import sys, os
sys.path.append(os.path.expanduser('../../../'))
sys.path.append(os.path.expanduser('../../../PyHEADTAIL/'))


import numpy as np
import pylab as pl
import mystyle as ms

n_segments = 5
B_multip = [0.]

# define machine for PyHEADTAIL
from PyHEADTAIL.particles.slicing import UniformBinSlicer
from machines_for_testing  import shortSPS
machine = shortSPS(n_segments = n_segments, machine_configuration = 'Q20-injection-like')

# remove synchrotron motion
machine.one_turn_map.remove(machine.longitudinal_map)


# compute sigma x and y
epsn_x = 2.5e-6
epsn_y = 2.5e-6

inj_optics = machine.transverse_map.get_injection_optics()
sigma_x = np.sqrt(inj_optics['beta_x'] * epsn_x / machine.betagamma)
sigma_y = np.sqrt(inj_optics['beta_y'] * epsn_y / machine.betagamma)


# generate a bunch
bunch = machine.generate_6D_Gaussian_bunch(n_macroparticles=300000, intensity=1.15e11, epsn_x=epsn_x, epsn_y=epsn_y, sigma_z=0.2)


# define apertures and Dh_sc to simulate headtail conditions
x_aper = 20 * sigma_x
y_aper = 20 * sigma_y
Dh_sc = 2 * x_aper / 128

# ecloud
import PyECLOUD.PyEC4PyHT as PyEC4PyHT
from PyHEADTAIL.particles.slicing import UniformBinSlicer
slicer = UniformBinSlicer(n_slices = 64, n_sigma_z = 3.)


init_unif_edens_flag = 1
init_unif_edens = 1e11
N_MP_ele_init = 100000
N_mp_max = N_MP_ele_init * 4.

nel_mp_ref_0 = init_unif_edens * 4 * x_aper * y_aper / N_MP_ele_init


ecloud = PyEC4PyHT.Ecloud(L_ecloud=machine.circumference / machine.transverse_map.n_segments, slicer=slicer,
				Dt_ref=25e-12, pyecl_input_folder='./drift_sim',
				x_aper=x_aper, y_aper=y_aper, Dh_sc=Dh_sc,
				init_unif_edens_flag=init_unif_edens_flag,
				init_unif_edens=init_unif_edens,
				N_mp_max=N_mp_max,
				nel_mp_ref_0=nel_mp_ref_0,
				B_multip=B_multip)
machine.install_after_each_transverse_segment(ecloud)


# set number of turns
n_turns = 128
n_record = 1500

# prepare storage for particles cohordinates
x_i = np.empty((n_record, n_turns))
xp_i = np.empty((n_record, n_turns))
y_i = np.empty((n_record, n_turns))
yp_i = np.empty((n_record, n_turns))

# track and store
for i in range(n_turns):
    machine.track(bunch)#, verbose=True)

    print 'Turn', i
    sys.stdout.flush()

    x_i[:,i] = bunch.x[:n_record]
    xp_i[:,i] = bunch.xp[:n_record]
    y_i[:,i] = bunch.y[:n_record]
    yp_i[:,i] = bunch.yp[:n_record]
print '\nDONE'

from tune_analysis import tune_analysis
qx_i, qy_i, qx_centroid, qy_centroid = tune_analysis(x_i, xp_i, y_i, yp_i)

pl.close('all')
pl.figure(1)
pl.plot(np.mean(x_i, axis=0))
pl.plot(np.mean(y_i, axis=0))


pl.figure(2)
pl.plot(np.abs(qx_i), np.abs(qy_i), '.')
pl.plot([np.modf(machine.Q_x)[0]], [np.modf(machine.Q_y)[0]], 'r.')
#plt.plot([np.abs(qx_centroid)], [np.abs(qy_centroid)], 'og')
pl.xlabel('$Q_x$');pl.ylabel('$Q_y$')
pl.axis('equal')

pl.show()

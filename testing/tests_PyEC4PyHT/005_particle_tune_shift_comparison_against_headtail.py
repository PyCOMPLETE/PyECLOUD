import sys, os
BIN=os.path.expanduser('../')
sys.path.append(BIN)


import numpy as np
import pylab as pl
import mystyle as ms


#~ filename = 'headtail_for_test/test_protons/shortSPS_Q20_proton_check_tune_spread_prb.dat'
#~ N_kicks = 5
#~ n_segments = 5
#~ B_multip = [0.]


filename = 'headtail_for_test/test_protons/shortSPS_Q20_proton_check_tune_spread_dipole_prb.dat'
N_kicks = 5
n_segments = 5
B_multip = [0.5]




n_record = 1000

n_part_per_turn = 5000

print 'Loading HEADTAIL data...'
appo = np.loadtxt(filename)

parid = np.reshape(appo[:,0], (-1, n_part_per_turn))[::N_kicks,:]
x = np.reshape(appo[:,1], (-1, n_part_per_turn))[::N_kicks,:]
xp = np.reshape(appo[:,2], (-1, n_part_per_turn))[::N_kicks,:]
y = np.reshape(appo[:,3], (-1, n_part_per_turn))[::N_kicks,:]
yp =np.reshape(appo[:,4], (-1, n_part_per_turn))[::N_kicks,:]
z = np.reshape(appo[:,5], (-1, n_part_per_turn))[::N_kicks,:]
zp = np.reshape(appo[:,6], (-1, n_part_per_turn))[::N_kicks,:]
n_turns = len(x[:,0])

print 'Done!'



# define machine for PyHEADTAIL
from PyHEADTAIL.particles.slicing import UniformBinSlicer
from shortSPS import shortSPS
machine = shortSPS(n_segments = n_segments, machine_configuration = 'Q20-injection-like')

# remove synchrotron motion
machine.one_turn_map.remove(machine.longitudinal_map)

# compute tunes from headtail data
print 'Tune analysis for headtail...' 
from tune_analysis import tune_analysis
qx_ht, qy_ht, qx_centroid_ht, qy_centroid_ht  = tune_analysis(None, machine, x[:, :n_record].T,
			xp[:, :n_record].T, y[:, :n_record].T, yp[:, :n_record].T)


# compute sigma x and y
epsn_x = 2.5e-6
epsn_y = 2.5e-6

sigma_x = np.sqrt(machine.beta_x[0]*epsn_x/machine.betagamma)
sigma_y = np.sqrt(machine.beta_y[0]*epsn_y/machine.betagamma)


				
# generate a bunch 
bunch = machine.generate_6D_Gaussian_bunch(n_macroparticles=300000, intensity=1.15e11, epsn_x=epsn_x, epsn_y=epsn_y, sigma_z=0.2)

# replace first particles with HEADTAIL ones
bunch.x[:n_part_per_turn] = x[0,:]
bunch.xp[:n_part_per_turn] = xp[0,:]
bunch.y[:n_part_per_turn] = y[0,:]
bunch.yp[:n_part_per_turn] = yp[0,:]
bunch.z[:n_part_per_turn] = z[0,:]
bunch.dp[:n_part_per_turn] =zp[0,:]


# define apertures and Dh_sc to simulate headtail conditions
x_aper  = 20*sigma_x
y_aper  = 20*sigma_y
Dh_sc = 2*x_aper/128/2.

# ecloud
import PyECLOUD.PyEC4PyHT as PyEC4PyHT
from PyHEADTAIL.particles.slicing import UniformBinSlicer
slicer = UniformBinSlicer(n_slices = 64, n_sigma_z = 3.)


init_unif_edens_flag=1
init_unif_edens=1e11
N_MP_ele_init = 100000
N_mp_max = N_MP_ele_init*4.

nel_mp_ref_0 = init_unif_edens*4*x_aper*y_aper/N_MP_ele_init


ecloud = PyEC4PyHT.Ecloud(L_ecloud=machine.circumference/machine.n_segments, slicer=slicer, 
				Dt_ref=25e-12, pyecl_input_folder='./drift_sim',
				x_aper=x_aper, y_aper=y_aper, Dh_sc=Dh_sc,
				init_unif_edens_flag=init_unif_edens_flag,
				init_unif_edens=init_unif_edens, 
				N_MP_ele_init=N_MP_ele_init,
				N_mp_max=N_mp_max,
				nel_mp_ref_0=nel_mp_ref_0,
				B_multip=B_multip)
machine.install_after_each_transverse_segment(ecloud)




# prepare storage for particles cohordinates
x_i = np.empty((n_record, n_turns))
xp_i = np.empty((n_record, n_turns))
y_i = np.empty((n_record, n_turns))
yp_i = np.empty((n_record, n_turns))

# track and store
for i in range(n_turns):    
    machine.track(bunch, verbose=True)
    
    print 'Turn', i
    sys.stdout.flush()
    
    x_i[:,i] = bunch.x[:n_record]
    xp_i[:,i] = bunch.xp[:n_record]
    y_i[:,i] = bunch.y[:n_record]
    yp_i[:,i] = bunch.yp[:n_record]
print '\nDONE'

from tune_analysis import tune_analysis
qx_i, qy_i, qx_centroid, qy_centroid  = tune_analysis(bunch, machine, x_i, xp_i, y_i, yp_i)

pl.close('all')
ms.mystyle(fontsz=14)
pl.figure(1)
sp1 = pl.subplot(2,1,1)
pl.plot(np.mean(x_i, axis=0), '.-b', markersize=5, linewidth=2, label='PyHT')
pl.plot(np.mean(x, axis=1),'.-r', markersize=5, linewidth=2, label='HT')
pl.ylabel('<x>')
pl.grid('on')
ms.sciy()
pl.legend(prop={'size':14})
pl.subplot(2,1,2, sharex=sp1)
pl.plot(np.mean(y_i, axis=0), '.-b', markersize=5, linewidth=2, label='PyHT')
pl.plot(np.mean(y, axis=1),'.-r', markersize=5, linewidth=2, label='HT')
pl.xlabel('Turn'); pl.ylabel('<y>')
pl.grid('on')
ms.sciy()
pl.savefig(filename.split('_prb.dat')[0]+'_centroids.png', dpi=200)


pl.figure(2)
pl.plot(np.abs(qx_i), np.abs(qy_i), '.', label='PyHT', markersize=3)
pl.plot(np.abs(qx_ht), np.abs(qy_ht), '.r', label='HT', markersize=3)
pl.plot([np.modf(machine.Q_x)[0]], [np.modf(machine.Q_y)[0]], 'go')
pl.xlabel('$Q_x$');pl.ylabel('$Q_y$')
pl.legend(prop={'size':14})
pl.grid('on')
pl.axis('equal')
pl.savefig(filename.split('_prb.dat')[0]+'_footprint.png', dpi=200)


pl.figure(3)
pl.subplot(2,1,1)
pl.plot(bunch.z[:n_record],np.abs(qx_i)-np.modf(machine.Q_x)[0], '.', markersize=3, label='PyHT')
pl.plot(z[-1,:n_record].T,np.abs(qx_ht)-np.modf(machine.Q_x)[0], '.r', markersize=3, label='HT')
pl.ylabel('$\Delta Q_x$')
pl.grid('on')
pl.legend(prop={'size':14})
pl.subplot(2,1,2)
pl.plot(bunch.z[:n_record],np.abs(qy_i)-np.modf(machine.Q_x)[0], '.', markersize=3)
pl.plot(z[-1,:n_record].T,np.abs(qy_ht)-np.modf(machine.Q_x)[0], '.r', markersize=3)
pl.ylabel('$\Delta Q_y$')
pl.xlabel('z [m]')
pl.grid('on')
pl.savefig(filename.split('_prb.dat')[0]+'_Q_vs_z.png', dpi=200)
pl.show()

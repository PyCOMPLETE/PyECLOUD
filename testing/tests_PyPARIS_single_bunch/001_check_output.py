import sys
sys.path.append('../../')

import numpy as np
import pylab as pl

import PyPARIS.myfilemanager as mfm
import mystyle as ms

import time

n_segments = 3
n_part_per_turn = 5000
N_turns = 8

filename = '../tests_PyEC4PyHT/headtail_for_test/test_protons/SPS_Q20_proton_check_dipole_3kicks_20150212_prb.dat'
appo = np.loadtxt(filename)

parid = np.reshape(appo[:,0], (-1, n_part_per_turn))[::n_segments,:]
x = np.reshape(appo[:,1], (-1, n_part_per_turn))[::n_segments,:]
xp = np.reshape(appo[:,2], (-1, n_part_per_turn))[::n_segments,:]
y = np.reshape(appo[:,3], (-1, n_part_per_turn))[::n_segments,:]
yp =np.reshape(appo[:,4], (-1, n_part_per_turn))[::n_segments,:]
z = np.reshape(appo[:,5], (-1, n_part_per_turn))[::n_segments,:]
zp = np.reshape(appo[:,6], (-1, n_part_per_turn))[::n_segments,:]

rms_err_x_list = []
rms_err_y_list = []
for i_turn in range(N_turns-1):
	
	ob = mfm.object_with_arrays_and_scalar_from_h5('particles_at_turn_%d.h5'%i_turn) 

	# sort id and momenta after track

	id_after = ob.id_after
	xp_after = ob.xp_after
	yp_after = ob.yp_after
	z_after = ob.z_after
	id_before = ob.id_before
	xp_before = ob.xp_before
	yp_before = ob.yp_before

	fig = pl.figure(100, figsize=(18,10))
	pl.clf()
	sp1 = pl.subplot(2,3,1)
	pl.plot(z_after,  (100*np.abs((xp_after-xp_before)-(xp[i_turn+1, :]-xp[0, :]))/np.std(xp[i_turn+1, :]-xp[0, :])), '.r', markersize = 2)
	pl.grid('on')
	pl.ylabel('''error($\Delta$x')/$\Delta$x' [%]''')
	pl.subplot(2,3,4, sharex=sp1)
	pl.plot(z_after,  (xp[i_turn+1, :]-xp[0, :]), '.r', label='HT', markersize = 2)
	pl.plot(z_after,  (xp_after-xp_before), '.b', label='PyHT', markersize = 2)
	ms.sciy()
	pl.grid('on')
	pl.ylabel("$\Delta$x'")
	pl.xlabel('z[m]')
	pl.legend(prop = {'size':14})
	pl.subplot(2,3,3)
	pl.hist((100*np.abs((xp_after-xp_before)-(xp[i_turn+1, :]-xp[0, :]))/np.std(xp[i_turn+1, :]-xp[0, :])), bins=100, range=(0,100))
	pl.xlabel('Error_x [%]')
	pl.ylabel('Occurrences')
	pl.xlim(0,100)
	rms_err_x = np.std(100*np.abs((xp_after-xp_before)-(xp[i_turn+1, :]-xp[0, :]))/np.std(xp[i_turn+1, :]-xp[0, :]))

	sp1 = pl.subplot(2,3,2)
	pl.plot(z_after,  (100*np.abs((yp_after-yp_before)-(yp[i_turn+1, :]-yp[0, :]))/np.std(yp[i_turn+1, :]-yp[0, :])), '.r', markersize = 2)
	pl.grid('on')
	pl.ylabel('''error($\Delta$y')/$\Delta$y' [%]''')
	pl.subplot(2,3,5, sharex=sp1)
	pl.plot(z_after,  (yp[i_turn+1, :]-yp[0, :]), '.r', label='HT', markersize = 2)
	pl.plot(z_after,  (yp_after-yp_before), '.b', label='PyHT', markersize = 2)
	ms.sciy()
	pl.grid('on')
	pl.ylabel("$\Delta$y'")
	pl.xlabel('z[m]')
	pl.legend(prop = {'size':14})
	pl.subplot(2,3,6)
	pl.hist((100*np.abs((yp_after-yp_before)-(yp[i_turn+1, :]-yp[0, :]))/np.std(yp[i_turn+1, :]-yp[0, :])), bins=100, range=(0,100))
	pl.xlim(0,100)
	pl.xlabel('Error_y [%]')
	pl.ylabel('Occurrences')
	rms_err_y = np.std(100*np.abs((yp_after-yp_before)-(yp[i_turn+1, :]-yp[0, :]))/np.std(yp[i_turn+1, :]-yp[0, :]))


	pl.suptitle('Turn %d rms_err_x = %e rms_err_y = %e'%(i_turn, rms_err_x, rms_err_y))

	pl.savefig(filename.split('_prb.dat')[0]+'_%02d.png'%i_turn, dpi=150)

	rms_err_x_list.append(rms_err_x)
	rms_err_y_list.append(rms_err_y)
	
	pl.ion()
	pl.draw()
	time.sleep(1.)
	
pl.figure(1000)
pl.plot(rms_err_x_list, '.-', markersize = 10, linewidth=2, label='x')
pl.plot(rms_err_y_list, '.-', markersize = 10, linewidth=2, label='y')
pl.grid('on')
pl.ylabel('''Relative r.m.s. error [%]''')
pl.xlabel('Turn')
pl.legend()
pl.savefig(filename.split('_prb.dat')[0]+'_errors.png', dpi=200)
pl.show()

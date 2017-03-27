import sys, os
sys.path.append(os.path.expanduser('../../../'))
sys.path.append(os.path.expanduser('../../../PyHEADTAIL/'))

print ''
print 'This script generates h5 file with tunes from Headtail output files.'
print 'Skip this if h5 is already available.'
print ''

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
dp = np.reshape(appo[:,6], (-1, n_part_per_turn))[::N_kicks,:]
n_turns = len(x[:,0])

print 'Done!'

# compute tunes from headtail data
print 'Tune analysis for headtail...' 
from tune_analysis import tune_analysis
qx_ht, qy_ht, qx_centroid_ht, qy_centroid_ht  = tune_analysis(x[:, :n_record].T, xp[:, :n_record].T, y[:, :n_record].T, yp[:, :n_record].T)

x0_HT = x[0,:]
xp0_HT = xp[0,:]
y0_HT = y[0,:]
yp0_HT = yp[0,:]
z0_HT = z[0,:]
dp0_HT =dp[0,:]


# Save
import h5py
dict_save = {\
'qx_ht':qx_ht,
'qy_ht':qy_ht,
'qx_centroid_ht':qx_centroid_ht,
'qy_centroid_ht':qy_centroid_ht,
'x0_HT':x0_HT,
'xp0_HT':xp0_HT,
'y0_HT':y0_HT,
'yp0_HT':yp0_HT,
'z0_HT':z0_HT,
'dp0_HT':dp0_HT,
'n_turns':n_turns}

with h5py.File('footprint_HT.h5', 'w') as fid:
        for kk in dict_save.keys():
                fid[kk] = dict_save[kk]


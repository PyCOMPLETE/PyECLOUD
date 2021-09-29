import numpy as np

import PyECLOUD.myfilemanager as mfm

filepath = '../../testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns/Pyecltest_angle3D_ref.mat'
filepath = 'Pyecltest.mat' 

dd = mfm.myloadmat(filepath)

fields = sorted(dd.keys())
fsizes = []

for ff in fields:
    fsizes.append(dd[ff].nbytes)

ind_sorted = np.argsort(fsizes)

for ind in ind_sorted:
    print(('%s : %f kB'%(fields[ind], fsizes[ind]/1e3)))

import matplotlib.pyplot as plt
plt.figure(1)
plt.plot(dd['t'], dd['Nel_timep'])


plt.show()

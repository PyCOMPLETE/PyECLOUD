import PyECLOUD.myfilemanager as mfm
import numpy as np
ob = mfm.myloadmat_to_obj('Pyecltest.mat')

import matplotlib.pyplot as plt
plt.close('all')

fig1 = plt.figure(1)
plt.plot(ob.xg_hist, np.sum(ob.energ_eV_impact_hist, axis=0))

fig2 = plt.figure(2, figsize=(8,6*1.5))
sp1 = plt.subplot(3,1,1)
sp1.plot(ob.t, ob.Nel_timep)
sp2 = plt.subplot(3,1,2, sharex=sp1)
sp2.plot(ob.t, ob.Qpatch_ave)
sp3 = plt.subplot(3,1,3, sharex=sp1)
sp3.plot(ob.t, ob.lam_t_array)
plt.show()

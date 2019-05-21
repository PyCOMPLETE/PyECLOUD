import numpy as np
import matplotlib.pyplot as plt

import PyECLOUD.myfilemanager as mfm

i_det = 6

ob = mfm.myloadmat_to_obj('Pyecltest.mat')

plt.close('all')
fig1 = plt.figure(1)
sp1 = fig1.add_subplot(1,1,1)
sp1.plot(ob.t/ob.b_spac, ob.Nel_timep)

fig2 = plt.figure()
sp21 = fig2.add_subplot(2,1,1)
mask_det = np.logical_and(ob.t>=i_det*ob.b_spac, ob.t<(i_det+1)*ob.b_spac) 
t_det = ob.t[mask_det]
t_det -= t_det[0]

sp21.plot(t_det, ob.Nel_timep[mask_det]/ob.Nel_timep[mask_det][0])
plt.show()

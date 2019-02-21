import PyECLOUD.myfilemanager as mfm
import numpy as np
from scipy.constants import e as qe

ob = mfm.myloadmat_to_obj('Pyecltest.mat')

import matplotlib.pyplot as plt
plt.close('all')

fig2 = plt.figure(2, figsize=(8,6*1.5))
sp1 = plt.subplot(4,1,1)
sp1.plot(ob.t/ob.b_spac, ob.Nel_timep)
sp1.set_ylabel('Number of e- [m^-1]')
sp2 = plt.subplot(4,1,2, sharex=sp1)
sp2.plot(ob.t/ob.b_spac, ob.Qpatch_ave)
sp2.set_ylabel('Q on the patch [C/m^2]')
sp3 = plt.subplot(4,1,3, sharex=sp1)
sp3.plot(ob.t/ob.b_spac, ob.sey_at_emax_patch)
sp3.set_ylabel('SEY at Emax (patch)')
sp4 = plt.subplot(4,1,4, sharex=sp1)
sp4.plot(ob.t/ob.b_spac, ob.lam_t_array)

fig2.subplots_adjust(
    top=0.95,
    bottom=0.05,
    left=0.125,
    right=0.9,
    hspace=0.295,
    wspace=0.2)

fig1 = plt.figure(1, figsize=(8*2., 6))
spprof = plt.subplot(1,2,1)
spprof.plot(ob.xg_hist, np.sum(ob.energ_eV_impact_hist, axis=0))
spprof.set_ylabel('Heat load [a.u.]')
spprof.set_xlabel('x [m]')


hl_left = np.sum(ob.energ_eV_impact_hist[:, ob.xg_hist<0], axis=1)
hl_right = np.sum(ob.energ_eV_impact_hist[:, ob.xg_hist>0], axis=1)
splr = plt.subplot(1,2,2) 
splr.plot(hl_left, 'b.-')
splr.plot(hl_right, 'r.-')
splr.set_xlabel('Bunch passage')
splr.set_ylabel('Heat load [a.u.]')

# crosscheck current on patch
mask_patch = ob.flag_charging>0

nel_impact_on_patch = np.sum(ob.nel_hist_impact_seg[:, mask_patch], axis=1)
nel_emit_on_patch = np.sum(ob.nel_hist_emit_seg[:, mask_patch], axis=1)

patch_area = np.sum(ob.L_edg[mask_patch])

accumulated_charge_m2 = -qe*np.cumsum(nel_impact_on_patch - nel_emit_on_patch)/patch_area
sp2.plot(accumulated_charge_m2)

plt.show()

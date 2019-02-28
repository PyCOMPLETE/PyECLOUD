import PyECLOUD.myfilemanager as mfm
import PyECLOUD.mystyle as ms

import numpy as np
from scipy.constants import e as qe

fname = 'Pyecltest.mat'
ob = mfm.myloadmat_to_obj(fname)

compare_against_charge_in_chamber = False
plot_charge_from_post_processing = False

import matplotlib.pyplot as plt
plt.close('all')
ms.mystyle_arial(fontsz=16)

fig2 = plt.figure(2, figsize=(8,6*1.5))
fig2.set_facecolor('w')
sp1 = plt.subplot(3,1,1)
sp1.plot(ob.t/ob.b_spac, ob.Nel_timep)
sp1.set_ylabel('Number of e-\n[m^-1]')
sp2 = plt.subplot(3,1,2, sharex=sp1)
sp2.plot(ob.t/ob.b_spac, ob.Qpatch_ave)
sp2.set_ylabel('Q on the patch\n[C/m^2]')
sp3 = plt.subplot(3,1,3, sharex=sp1)
sp3.plot(ob.t/ob.b_spac, ob.sey_at_emax_patch)
sp3.set_ylabel('SEY at Emax\n(patch)')
# sp4 = plt.subplot(4,1,4, sharex=sp1)
# sp4.plot(ob.t/ob.b_spac, ob.lam_t_array)

from matplotlib.ticker import MaxNLocator
for sp in [sp1, sp2, sp3]:
    sp.yaxis.set_major_locator(MaxNLocator(5))
    sp.grid(True)
sp3.set_xlabel('Time/(25 ns)')

fig2.subplots_adjust(
    top=0.95,
    bottom=0.09,
    left=0.15,
    right=0.94,
    hspace=0.38,
    wspace=0.2)

fig1 = plt.figure(1, figsize=(8*2., 6))
fig1.set_facecolor('w')
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

fig1.suptitle(fname)
fig2.suptitle(fname)

# crosscheck current on patch
mask_patch = ob.flag_charging>0

nel_impact_on_patch = np.sum(ob.nel_hist_impact_seg[:, mask_patch], axis=1)
nel_emit_on_patch = np.sum(ob.nel_hist_emit_seg[:, mask_patch], axis=1)

patch_area = np.sum(ob.L_edg[mask_patch])

accumulated_charge_m2 = -qe*np.cumsum(nel_impact_on_patch - nel_emit_on_patch)/patch_area
if plot_charge_from_post_processing:
    sp2.plot(accumulated_charge_m2)


if compare_against_charge_in_chamber:
    set_patch_Vx = set(list(ob.Vx[1:][mask_patch]) + list(ob.Vx[:-1][mask_patch]))
    
    min_x_patch = np.min(list(set_patch_Vx))
    max_x_patch = np.max(list(set_patch_Vx))
    
    mask_xg_patch = np.logical_and(ob.xg_hist>min_x_patch, ob.xg_hist<max_x_patch)
    
    Q_chamb_patch = np.sum(ob.nel_hist[:, mask_xg_patch], axis=1)*qe/patch_area
    
    sp2.plot(Q_chamb_patch-Q_chamb_patch[0])
    
    sp1.plot(ob.t/ob.b_spac, ob.Qpatch_ave*patch_area/qe)

fig200 = plt.figure(200)
fig200.set_facecolor('w')
spchm = fig200.add_subplot(1,1,1)
spchm.plot(ob.Vx, ob.Vy, 'b.-', linewidth=2, markersize=10)
for ii in np.where(mask_patch)[0]:
    spchm.plot(ob.Vx[ii:ii+2], ob.Vy[ii:ii+2], 'r.-',
            linewidth=2, markersize=10)
spchm.axis('equal')


plt.show()

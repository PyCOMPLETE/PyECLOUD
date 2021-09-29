import numpy as np
import matplotlib.pyplot as plt

import PyECLOUD.myfilemanager as mfm
import PyECLOUD.mystyle as ms
obfol = mfm.myloadmat_to_obj('followed_electrons.mat')
ob = mfm.myloadmat_to_obj('Pyecltest.mat')

i_obs = 1

plt.close('all')
ms.mystyle(fontsz=14, traditional_look=False)

Dz_all = []
nel_all = []

nele_no_nan = obfol.nel.copy()
nele_no_nan[np.isnan(obfol.nel)] = 0.

for i_ele in range(obfol.nel.shape[1]):

    nel_ele = nele_no_nan[:, i_ele]
    z_ele = obfol.z[:, i_ele]

    i_changes = np.where(np.abs(np.diff(nel_ele)) > 0.)[0]

    for ii in range(1, len(i_changes)-1):
        z_part = z_ele[i_changes[ii]+1: i_changes[ii+1]-1]
        if len(z_part) > 0 and not np.isnan(z_part[0]):
            Dz_all.append(np.max(z_part)-np.min(z_part))
            nel_all.append(nel_ele[i_changes[ii]+1])

    if i_ele == i_obs:
        fig1 = plt.figure(1)
        ax1 = fig1.add_subplot(2,1,1)
        ax2 = fig1.add_subplot(2,1,2, sharex=ax1)

        ax1.plot(nel_ele)
        ax1.plot(i_changes, nel_ele[i_changes], '.')
        ax2.plot(z_ele)

hist, bin_edges = np.histogram(Dz_all, bins=100, range=(0, 0.4), weights=nel_all)

fig2 = plt.figure(2)
ax21 = fig2.add_subplot(1,1,1)
#ax21.plot(0.5*(bin_edges[1:]+bin_edges[:-1]), hist)
ax21.bar(x=0.5*(bin_edges[1:]+bin_edges[:-1]), height=hist,
        align='center', width=np.mean(np.diff(bin_edges)))
ax21.set_xlim(left=0.)
ax21.set_ylim(bottom=0.)
ax21.set_ylabel('Occurrencies')
ax21.set_xlabel('Longitudinal displacement [m]')
ax21.grid(True, linestyle=':', alpha=.5)
plt.show()

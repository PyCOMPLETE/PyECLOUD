import numpy as np
import matplotlib.pyplot as plt

import PyECLOUD.myfilemanager as mfm

ob = mfm.myloadmat_to_obj('followed_electrons.mat')

i_obs = 1

Dz_all = []

nele_no_nan = ob.nel.copy()
nele_no_nan[np.isnan(ob.nel)] = 0.

for i_ele in range(ob.nel.shape[1]):

    nel_ele = nele_no_nan[:, i_ele]
    z_ele = ob.z[:, i_ele]

    i_changes = np.where(np.abs(np.diff(nel_ele)) > 0.)[0]

    for ii in range(1, len(i_changes)-1):
        z_part = z_ele[i_changes[ii]+1: i_changes[ii+1]-1]
        if len(z_part) > 0 and not np.isnan(z_part[0]):
            Dz_all.append(np.max(z_part)-np.min(z_part))

    if i_ele == i_obs:
        plt.close('all')
        fig1 = plt.figure(1)
        ax1 = fig1.add_subplot(2,1,1)
        ax2 = fig1.add_subplot(2,1,2, sharex=ax1)

        ax1.plot(nel_ele)
        ax1.plot(i_changes, nel_ele[i_changes], '.')
        ax2.plot(z_ele)

plt.show()

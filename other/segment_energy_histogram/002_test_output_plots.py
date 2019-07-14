
import matplotlib.pyplot as plt
pl = plt
import numpy as np

import PyECLOUD.myloadmat_to_obj as mlo
import PyECLOUD.mystyle as ms

pl.close('all')
ms.mystyle_arial(fontsz=16, dist_tick_lab=5)
dpiset = 200

ob = mlo.myloadmat_to_obj('./Pyecltest.mat')


mid_x = 0.5*(ob.Vx[1:]+ob.Vx[:-1])

Dx_plot = 3e-3
x_spectrum = np.arange(0., np.max(mid_x), Dx_plot)
spectra = []
for ii, xx in enumerate(x_spectrum):
    mask = np.abs(mid_x - xx) < Dx_plot # It is meant to overlap a bit
    spec = np.sum(np.sum(ob.En_hist_seg[mask, :, :], axis=0), axis=0)
    spec /= np.sum(spec)
    spectra.append(spec)

fig10 = plt.figure(10)
fig10.set_facecolor('w')
ax10 = fig10.add_subplot(1,1,1)
for ii, xx in enumerate(x_spectrum):
    col = ms.colorprog(ii, len(x_spectrum))
    ax10.semilogy(ob.En_g_hist, spectra[ii], 
            color=col, linewidth=2., label='%.1f mm'%(xx*1e3))

fig11 = plt.figure(11)
ax11 = fig11.add_subplot(1,1,1)
ax11.semilogy(ob.xg_hist*1e3, np.sum(ob.nel_impact_hist_tot, axis=0),
        linewidth=2.)

ax10.legend(prop={'size':16})
ax10.set_ylim(bottom=1e-4)
ax10.grid(True)
ax10.set_ylabel('Normalized energy spectrum')
ax10.set_xlabel('Energy [eV]')

fig10.subplots_adjust(bottom=.14)

fig = pl.figure(1)
pl.plot(ob.t, ob.Nel_timep, linewidth=2)
pl.xlabel('Time [s]')
pl.ylabel('Number of $e^-$ per unit length [$m^{-1}$]')
ms.scix()
pl.grid('on')
pl.suptitle('Var. name: Nel_timep\nNumber of electrons in the chamber at each time step')
pl.subplots_adjust(top=.82, bottom=.14)




pl.show()

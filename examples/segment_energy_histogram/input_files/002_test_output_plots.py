
import matplotlib.pyplot as plt
pl = plt
import numpy as np

import PyECLOUD.myloadmat_to_obj as mlo
import PyECLOUD.mystyle as ms

pl.close('all')
ms.mystyle_arial(fontsz=16)
dpiset = 200

ob = mlo.myloadmat_to_obj('./Pyecltest.mat')


mid_x = 0.5*(ob.Vx[1:]+ob.Vx[:-1])

Dx_plot = 2e-3
x_spectrum = np.arange(0., np.max(mid_x), Dx_plot)
spectra = []
for ii, xx in enumerate(x_spectrum):
    mask = np.abs(mid_x - xx) < Dx_plot # It is meant to overlap a bit
    spec = np.sum(np.sum(ob.En_hist_seg[mask, :, :], axis=0), axis=0)
    spec /= np.sum(spec)
    spectra.append(spec)

fig10 = plt.figure(10)
for ii, xx in enumerate(x_spectrum):
    col = ms.colorprog(ii, len(x_spectrum))
    plt.semilogy(ob.En_g_hist, spectra[ii], 
            color=col, linewidth=2., label='%.1f mm'%(xx*1e3))
plt.legend()

fig = pl.figure(1)
pl.plot(ob.t, ob.Nel_timep, linewidth=2)
pl.xlabel('Time [s]')
pl.ylabel('Number of $e^-$ per unit length [$m^{-1}$]')
ms.scix()
pl.grid('on')
pl.suptitle('Var. name: Nel_timep\nNumber of electrons in the chamber at each time step')
pl.subplots_adjust(top=.82, bottom=.14)




pl.show()

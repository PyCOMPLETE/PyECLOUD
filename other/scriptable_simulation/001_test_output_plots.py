# generate example figures with main PyECLOUD outputs
# Thanks for E. Belli, L. Mether and A. Romano for preparing this script

import sys
import os
BIN = os.path.expanduser("../../../../")
sys.path.append(BIN)

import pylab as pl
import numpy as np

import PyECLOUD.myloadmat_to_obj as mlo
import PyECLOUD.mystyle as ms

pl.close('all')
ms.mystyle_arial(fontsz=16)
dpiset = 200

ob = mlo.myloadmat_to_obj('./Pyecltest.mat')


#################################
# Variables saved per time step #
#################################

fig = pl.figure()
pl.plot(ob.t, ob.Nel_timep, linewidth=2)
pl.plot(ob.t, ob.Nelectrons, 'r--', linewidth=2)

pl.xlabel('Time [s]')
pl.ylabel('Number of $e^-$ per unit length [$m^{-1}$]')
ms.scix()
pl.grid('on')
pl.suptitle('Var. name: Nel_timep\nNumber of electrons in the chamber at each time step')
pl.subplots_adjust(top=.82, bottom=.14)

fig = pl.figure()
sp1 = pl.subplot(4,1,1)
sp1.plot(ob.t, ob.y_first_ele, '.-')
sp1.set_ylabel('y [m]')
sp3 = pl.subplot(4,1,2, sharex=sp1)
sp3.plot(ob.t, ob.relative_charge_first_ele, '.-')
sp3.set_ylabel('Nel/nel_ref')
sp4 = pl.subplot(4,1,3, sharex=sp1)
sp4.plot(ob.t, ob.ref_size, '.-')
sp4.set_ylabel('Ref MP size')
sp2 = pl.subplot(4,1,4, sharex=sp1)
sp2.plot(ob.t, ob.lam_t_array, '.-')
sp2.set_ylabel('Beam profile')
pl.show()

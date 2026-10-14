import sys, os
BIN = os.path.expanduser("../../../../")
sys.path.append(BIN)

import argparse
import matplotlib.pyplot as plt
import numpy as np
#from colorsys import hsv_to_rgb
import os
import PyECLOUD.myloadmat_to_obj as mlm
import matplotlib.gridspec as gridspec
import PyECLOUD.mystyle as ms

plt.close('all')
ms.mystyle_arial(fontsz=16)

ob = mlm.myloadmat_to_obj('./test_saving__iter0.mat')   # load dictionary of the current simulation

N_pass = ob.nel_hist.shape[0]

fig1 = plt.figure(1)
fig1.set_facecolor('w')
plt.pcolormesh(ob.xg_hist, ob.t_hist * 1e6, ob.nel_hist, shading='Gouraud')
plt.axis('tight')

x_bun = 0*ob.t_hist
x_bun[30:] += 3e-3
plt.plot(x_bun, ob.t_hist * 1e6, '.w', markersize=10.)
plt.xlabel('x [m]')
plt.ylabel('Time [us]')
fig1.subplots_adjust(bottom=.12)

plt.show()

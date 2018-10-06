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


ms.mystyle_arial(fontsz=14)

ob = mlm.myloadmat_to_obj('./test_saving__iter0.mat')   # load dictionary of the current simulation

N_pass = ob.nel_hist.shape[0]

plt.figure(1)
plt.pcolormesh(ob.xg_hist, range(N_pass), ob.nel_hist, shading='Gouraud')
plt.show()

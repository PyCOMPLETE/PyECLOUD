import sys, os
BIN = os.path.expanduser("../../../")
sys.path.append(BIN)

import argparse
import pylab as pl
import numpy as np
#from colorsys import hsv_to_rgb
import os
import PyECLOUD.myloadmat_to_obj as mlm
import matplotlib.gridspec as gridspec
import PyECLOUD.mystyle as ms

pl.close('all')

sim_folder = './'



myfontsz = 10
ms.mystyle_arial(fontsz=myfontsz)

dict_ref = mlm.myloadmat('../tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_multigrid/'+'/Pyecltest_angle3D_ref.mat') # load dictionary of the reference simulation
dict_curr = mlm.myloadmat('./test_saving.mat')   # load dictionary of the current simulation


out_var_ref = dict_ref.keys()       # returns the list of keys
out_var_curr = dict_curr.keys()

out_var_ref.sort()         # sort the keys in alphabetical order
out_var_curr.sort()

n_pass = 5  # reference passage

pl.close('all')

pl.figure(1000)
pl.plot(dict_ref['t'], dict_ref['Nel_timep'], '.-')
pl.plot(dict_curr['t']-10e-9, dict_curr['Nel_timep'], '.r-')

pl.figure(1001)
pl.semilogy(dict_ref['t'], dict_ref['En_kin_eV_time'], '.-')
pl.semilogy(dict_curr['t']-10e-9, dict_curr['En_kin_eV_time'], '.r-')

pl.figure(1002)
pl.plot(dict_ref['t'][:-1], np.diff(dict_ref['t'])/1e-12, '.-')
pl.plot(dict_curr['t'][:-1]-10e-9, np.diff(dict_curr['t'])/1e-12, '.r-')

pl.figure(1003)
pl.plot(dict_ref['N_mp_pass'], '.-')
pl.plot(dict_curr['N_mp_pass'], '.r-')

pl.figure(1004)
i_pass = 20
pl.plot(dict_ref['nel_hist'][i_pass, :], '.-')
pl.plot(dict_curr['nel_hist'][i_pass, :], '.r-')

pl.figure(1005)
pl.plot(dict_ref['En_hist'][i_pass, :], '.-')
pl.plot(dict_curr['En_hist'][i_pass, :], '.r-')

pl.figure(1006)
pl.plot(dict_ref['t'][:-1], dict_ref['En_imp_eV_time'][:-1]/np.diff(dict_ref['t']), '.-')
pl.plot(dict_curr['t'][:-1]-10e-9, dict_curr['En_imp_eV_time'][:-1]/np.diff(dict_curr['t']), '.r-')


pl.figure(1007)
pl.plot(dict_ref['energ_eV_impact_hist'][i_pass, :], '.-')
pl.plot(dict_curr['energ_eV_impact_hist'][i_pass, :], '.r-')

pl.figure(1008)
pl.semilogy(dict_ref['t'], dict_ref['lam_t_array'], '.-')
pl.semilogy(dict_curr['t']-10e-9, dict_curr['lam_t_array'], '.r-')


pl.show()
pl.show()




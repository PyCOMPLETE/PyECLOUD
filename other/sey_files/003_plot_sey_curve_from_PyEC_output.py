import numpy as np
import matplotlib.pyplot as pl
import sys
sys.path.append('../../../')

import PyECLOUD.myloadmat_to_obj as mlo
import PyECLOUD.mystyle as ms

pyecl_output_file = '../../testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns/Pyecltest_angle3D.mat'
#~ pyecl_output_file = '../../testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_change_s_and_E0/Pyecltest_angle3D.mat'
pyecl_output_file = '../../testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_seyfromfile/Pyecltest_angle3D.mat'


ob = mlo.myloadmat_to_obj(pyecl_output_file)

pl.close('all')
ms.mystyle_arial()

fig1 = pl.figure(1, figsize=(2 * 8, 6))
fig1.set_facecolor('w')
sp1 = fig1.add_subplot(1, 2, 1)
sp2 = fig1.add_subplot(1, 2, 2, sharex=sp1)
for i_ct, ct in enumerate(ob.sey_test_cos_theta):
    thiscol = ms.colorprog(i_ct, len(ob.sey_test_cos_theta))
    label = 'costheta=%.2f'%ct
    sp1.plot(ob.sey_test_E_impact_eV, ob.sey_test_del_true_mat[i_ct, :], color=thiscol, label=label)
    sp2.plot(ob.sey_test_E_impact_eV, ob.sey_test_del_elast_mat[i_ct, :], color=thiscol, label=label)

sp2.legend(loc='best', prop={'size': 14})
sp1.set_ylabel('Delta true')
sp2.set_ylabel('Delta elast')

for sp in [sp1, sp2]:
    sp.grid('on')
    sp.set_xlabel('Electron energy [eV]')


pl.show()


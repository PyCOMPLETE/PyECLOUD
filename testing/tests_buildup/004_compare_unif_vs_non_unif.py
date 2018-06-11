import os, sys

BIN = os.path.expanduser("../../../") #folder containing PyECLOUD, PyPIC, PyKLU
if BIN not in sys.path:
    sys.path.append(BIN)

import PyECLOUD.myloadmat_to_obj as mlm

import matplotlib.pyplot as plt

# obu = mlm.myloadmat_to_obj('LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns/Pyecltest_angle3D_ref.mat')
# obn = mlm.myloadmat_to_obj('LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_nonuniftime/Pyecltest_angle3D.mat')

obu = mlm.myloadmat_to_obj('LHC_ArcDipReal_450GeV_sey1.00_2.5e11ppb_bl_1.00ns_gas_ionization/Pyecltest_angle3D_ref.mat')
obn = mlm.myloadmat_to_obj('LHC_ArcDipReal_450GeV_sey1.00_2.5e11ppb_bl_1.00ns_gas_ionization/Pyecltest_angle3D.mat')

# obu = mlm.myloadmat_to_obj('LHC_ArcDipReal_6500GeV_sey_1.70_1.1e11ppb_b1_1.00ns/Pyecltest_angle3D_ref.mat')
# obn = mlm.myloadmat_to_obj('LHC_ArcDipReal_6500GeV_sey_1.70_1.1e11ppb_b1_1.00ns/Pyecltest_angle3D.mat')


plt.close('all')
plt.figure(1)
plt.plot(obu.t, obu.Nel_timep, '.-')
plt.plot(obn.t, obn.Nel_timep, '.r')
plt.show()

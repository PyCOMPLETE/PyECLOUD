import numpy as np
import matplotlib.pyplot as pl
import sys
sys.path.append('../../../')

from PyECLOUD.init import read_parameter_files, read_input_files_and_init_components
import PyECLOUD.mystyle as ms


pyecl_input_folder = '../../testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns/'
#~ pyecl_input_folder = '../../testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_change_s_and_E0/'
pyecl_input_folder = '../../testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_seyfromfile/'

beamtim,MP_e, dynamics,impact_man, pyeclsaver, \
    gas_ion_flag, resgasion, t_ion, \
    spacech_ele,t_sc_ON, photoem_flag, phemiss,\
    flag_presence_sec_beams, sec_beams_list, config_dict = read_input_files_and_init_components(
            pyecl_input_folder=pyecl_input_folder, skip_beam=True,
            skip_pyeclsaver=True, skip_spacech_ele=True)


n_rep = 10000
E_impact_eV_test = np.array(list(np.arange(0, 499., 1.)) + list(np.arange(500., 2000, 5)))
cos_theta_test = np.linspace(0, 1., 10)

del_true_mat, del_elast_mat = impact_man.extract_sey_curves(n_rep, E_impact_eV_test, cos_theta_test)

pl.close('all')
ms.mystyle_arial()

fig1 = pl.figure(1, figsize=(2*8,6))
fig1.set_facecolor('w')
sp1 = fig1.add_subplot(1,2,1)
sp2 = fig1.add_subplot(1,2,2, sharex=sp1)
for i_ct, ct in enumerate(cos_theta_test):
    thiscol = ms.colorprog(i_ct, len(cos_theta_test))
    label = 'costheta=%.2f'%ct
    sp1.plot(E_impact_eV_test, del_true_mat[i_ct, :], color=thiscol, label=label)
    sp2.plot(E_impact_eV_test, del_elast_mat[i_ct, :], color=thiscol, label=label)

sp2.legend(loc='best', prop={'size':14})
sp1.set_ylabel('Delta true')
sp2.set_ylabel('Delta elast')

for sp in [sp1, sp2]:
    sp.grid('on')
    sp.set_xlabel('Electron energy [eV]')


pl.show()


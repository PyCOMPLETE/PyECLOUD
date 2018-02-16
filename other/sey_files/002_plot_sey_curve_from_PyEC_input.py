import numpy as np
import matplotlib.pyplot as pl
import sys
sys.path.append('../../../') 

from PyECLOUD.init import read_parameter_files, read_input_files_and_init_components
import PyECLOUD.mystyle as ms


pyecl_input_folder = '../../testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_change_s_and_E0/'

#~ config_dict = read_parameter_files(pyecl_input_folder, skip_beam_files=True)

beamtim,MP_e, dynamics,impact_man, pyeclsaver, \
    gas_ion_flag, resgasion, t_ion, \
    spacech_ele,t_sc_ON, photoem_flag, phemiss,\
    flag_presence_sec_beams, sec_beams_list, config_dict = read_input_files_and_init_components(
            pyecl_input_folder=pyecl_input_folder, skip_beam=True,
            skip_pyeclsaver=True, skip_spacech_ele=True)
            
sey_mod = impact_man.sey_mod

E_impact_eV = np.array(list(np.arange(0, 499., 1.)) + list(np.arange(500., 2000, 5)))
cos_theta = np.linspace(0, 1., 10)
nel_impact = 1. + 0.*E_impact_eV

n_rep = 1000

del_true_mat = np.zeros((len(cos_theta), len(E_impact_eV)))
del_elast_mat = np.zeros((len(cos_theta), len(E_impact_eV)))

for i_ct, ct in enumerate(cos_theta):
    print('%d/%d'%(i_ct+1, len(cos_theta))) 
    for i_ene, Ene in enumerate(E_impact_eV):
        
        nel_emit, flag_elast, flag_truesec = sey_mod.SEY_process(nel_impact=np.ones(n_rep), 
                        E_impact_eV=Ene*np.ones(n_rep), costheta_impact=np.ones(n_rep)*ct, i_impact=np.array(n_rep*[0]))
                        
        del_true_mat[i_ct, i_ene] = np.mean(nel_emit)*float(np.sum(flag_truesec))/float(n_rep)
        del_elast_mat[i_ct, i_ene] = np.mean(nel_emit)*float(np.sum(flag_elast))/float(n_rep)
        

pl.close('all')
ms.mystyle_arial()

fig1 = pl.figure(1, figsize=(2*8,6))
fig1.set_facecolor('w')
sp1 = fig1.add_subplot(1,2,1)
sp2 = fig1.add_subplot(1,2,2, sharex=sp1)
for i_ct, ct in enumerate(cos_theta):
    thiscol = ms.colorprog(i_ct, len(cos_theta))
    label = 'costheta=%.2f'%ct
    sp1.plot(E_impact_eV, del_true_mat[i_ct, :], color=thiscol, label=label)
    sp2.plot(E_impact_eV, del_elast_mat[i_ct, :], color=thiscol, label=label)

sp2.legend(loc='best', prop={'size':14})  
sp1.set_ylabel('Delta true')
sp2.set_ylabel('Delta elast') 
        
for sp in [sp1, sp2]:
    sp.grid('on')
    sp.set_xlabel('Electron energy [eV]')
    

pl.show()


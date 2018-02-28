import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('../../../') 

import resample_sey as rss
import PyECLOUD.sec_emission_model_from_file as semf

fname_txt_file = 'sey_table_not_uniform.txt'
fname_unif_matfile = 'sey_unif.mat'
delta_E_unif_eV = 2.
range_for_extrap_eV = 300.

###  Loading table from lab measurements
data_table = np.loadtxt('sey_table_not_uniform.txt')
data_energy_eV = data_table[:, 0]
data_SEY = data_table[:, 1]

###  Resample data
dict_resampled = rss.resample_sey_data(energy_eV_samples=data_energy_eV, 
                    sey_true_samples=data_SEY, sey_elast_samples = data_SEY*0., 
                    uniform_dE=delta_E_unif_eV, range_extrapolate_right=range_for_extrap_eV)
                    
                    
### Save resampled data
import scipy.io as sio
sio.savemat(fname_unif_matfile,dict_resampled,oned_as='row')

### Test extrapolation using run-time SEY routine

# Create SecEmission object (as used in simulation)
se_obj = semf.SEY_model_from_file(sey_file=dict_resampled, flag_factor_costheta=True)

# Test impacts
N_ele_test = 10000
ene_eV_max_test = 3000.
ene_test_eV = np.linspace(0., ene_eV_max_test, N_ele_test)
nel_impact = ene_test_eV*0.+1.
costheta_impact = ene_test_eV*0.+1.

sey_test, refl, notrefl = se_obj.SEY_process(nel_impact=nel_impact,E_impact_eV=ene_test_eV, 
                            costheta_impact = costheta_impact, i_impact=None)
                            

# Plots!
plt.close('all')
plt.figure(1)
sp1 = plt.subplot(2,1,1) 
plt.plot(dict_resampled['energy_eV'], dict_resampled['sey_true'], 'r.-')
plt.plot(data_energy_eV, data_SEY, 'o')
plt.plot(ene_test_eV, sey_test, 'g')
sp2 = plt.subplot(2,1,2, sharex=sp1) 
plt.plot(dict_resampled['energy_eV'], dict_resampled['sey_elast'], 'r.-')
plt.plot(data_energy_eV, data_SEY*0, 'o')
#~ plt.plot(ene_test_eV, sey_test, 'g')


plt.show()



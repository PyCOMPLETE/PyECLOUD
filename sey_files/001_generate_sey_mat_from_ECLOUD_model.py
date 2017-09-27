import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('../') 

import sec_emission_model_ECLOUD as seme
import sec_emission_model_from_file as semf

Emax = 300.
del_max = 1.5
R0 = 0.7

costheta = 1.
E_samples_eV = np.linspace(0, 2000, 500)

dE_resample_eV = 1.
range_for_extrap_eV = 100

# Get curves from standard ECLOUD model
obec = seme.SEY_model_ECLOUD(Emax=Emax, del_max=del_max, R0=R0)
delta, ref_frac = seme.yield_fun2(E=E_samples_eV, costheta=costheta, 
          Emax=obec.Emax, del_max=obec.del_max, R0 = obec.R0, E0 = obec.E0, 
          s = obec.s)         
delta_true = delta*(1.-ref_frac)
delta_elast = delta*ref_frac


###  Resample data
import resample_sey as rss
dict_resampled = rss.resample_sey_data(energy_eV_samples=E_samples_eV, 
                    sey_true_samples=delta_true, sey_elast_samples = delta_elast, 
                    uniform_dE=dE_resample_eV, range_extrapolate_right=range_for_extrap_eV)
                    
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


plt.close('all')
plt.figure(1)
sp1 = plt.subplot(1,1,1) 
sp1.plot(E_samples_eV, delta, 'ob')
plt.plot(ene_test_eV, sey_test, 'b')
sp1.plot(E_samples_eV, delta_true, 'r')
sp1.plot(E_samples_eV, delta_elast, 'g')


plt.show()

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('../../../')

import PyECLOUD.sec_emission_model_ECLOUD as seme
import PyECLOUD.sec_emission_model_from_file as semf

Emax = 332.
del_max = 1.700000
R0 = 0.7

costheta = 1.

E_max_samples = 10000.#2000
N_samples = 100000 #500

dE_resample_eV = 0.2
range_for_extrap_eV = 100

E_samples_eV = np.linspace(0, E_max_samples, N_samples)


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

### Save resampled data
import scipy.io as sio
sio.savemat('sey_data_resampling_ecloud_model.mat',dict_resampled,oned_as='row')


# Create SecEmission object (as used in simulation)
se_obj = semf.SEY_model_from_file(sey_file=dict_resampled, flag_factor_costheta=True)


# Test impacts
N_ele_test = 10000
ene_eV_max_test = 10000.
ene_test_eV = np.linspace(0., ene_eV_max_test, N_ele_test)
nel_impact = ene_test_eV*0.+1.
costheta_impact = ene_test_eV*0.+1.

sey_test, refl, notrefl = se_obj.SEY_process(nel_impact=nel_impact,E_impact_eV=ene_test_eV,
                            costheta_impact = costheta_impact, i_impact=None)


E_test_refl_frac = list(np.linspace(0, 100, 100)) + list(np.linspace(100, 3000, 100))


N_this = 10000

delta_true_montecarlo = []
delta_elast_montecarlo = []
for E_this in E_test_refl_frac:
    nel_impact_this = np.ones(N_this)
    ene_test_eV_this = np.ones(N_this)*E_this
    sey_this, refl_this, notrefl_this = se_obj.SEY_process(nel_impact=nel_impact_this,E_impact_eV=ene_test_eV_this,
                                costheta_impact = costheta_impact, i_impact=None)
    N_refl_this = np.sum(refl_this)
    N_true_this = np.sum(notrefl_this)

    delta_true_this = sey_this[0] * float(N_true_this)/float(N_this)
    delta_elast_this = sey_this[0] * float(N_refl_this)/float(N_this)

    delta_true_montecarlo.append(delta_true_this)
    delta_elast_montecarlo.append(delta_elast_this)


plt.close('all')
plt.figure(1)
sp1 = plt.subplot(1,1,1)
sp1.plot(E_samples_eV, delta, 'ob', label='from ECLOUD mdl')
plt.plot(ene_test_eV, sey_test, 'b', label='from From-File mdl')
sp1.plot(E_samples_eV, delta_true, 'r', label='sample true')
sp1.plot(E_samples_eV, delta_elast, 'g', label='sample elast')

sp1.plot(E_test_refl_frac, delta_true_montecarlo, 'vr', label='MC true')
sp1.plot(E_test_refl_frac, delta_elast_montecarlo, 'vg', label='MC elast')

sp1.legend(loc='best')
sp1.grid('on')

plt.show()

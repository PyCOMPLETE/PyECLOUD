import sys
if '../..' not in sys.path:
    sys.path.append('../..')
import numpy as np
import matplotlib.pyplot as plt
import sec_emission_model_furman_pivi as fp
import mystyle as ms
from impact_management_class import impact_management
from geom_impact_ellip import ellip_cham_geom_object


plt.close('all')
ms.mystyle(35)
linewid = 2

me = 9.10938356e-31

furman_pivi_surface_tweak = {'conserve_energy': True,
                             'exclude_rediffused': True,
                             'choice': 'poisson',
                             'M_cut': 10,
                             'p_n': np.array([1.21963859, 1.66070543, 1.21935223, 1.09987752, 4.28158656, 1.02052557, 1.0247471, 1.02307995, 29.93491271, 1.02045612]),
                             'eps_n': np.array([7.44033631e+00, 2.47339424e+00, 7.45004962e+00, 1.63618903e+01, 4.97986255e-01, 7.96170380e+01, 6.60354258e+01, 7.08053955e+01, 5.64779654e-02, 7.98873331e+01]),
                             # Parameters for backscattered electrons
                             'p1EInf': 0.002158,  # Changed this
                             'p1Ehat': 0.709633,  # Changed this
                             'eEHat': 0.,
                             'w': 46.028959,  # Changed this
                             'p': 0.468907,  # Changed this
                             'e1': 0.,  # Changed this
                             'e2': 2.,
                             'sigmaE': 2.,
                             # Parameters for rediffused electrons
                             'p1RInf': 0.2,
                             'eR': 0.041,
                             'r': 0.104,
                             'q': 0.5,
                             'r1': 0.26,
                             'r2': 2.,
                             # Parameters for true secondaries
                             'deltaTSHat': 1.8848,
                             'eHat0': 332.,
                             's': 1.35,
                             't1': 0.706340,  # Changed this
                             't2': 0.715223,  # Changed this
                             't3': 0.7,
                             't4': 1.,
                             }

furman_pivi_surface_LHC = {'conserve_energy': True,
                           'exclude_rediffused': False,
                           'choice': 'poisson',
                           'M_cut': 10,
                           'p_n': np.array([2.5, 3.3, 2.5, 2.5, 2.8, 1.3, 1.5, 1.5, 1.5, 1.5]),
                           'eps_n': np.array([1.5, 1.75, 1., 3.75, 8.5, 11.5, 2.5, 3., 2.5, 3.]),
                           # Parameters for backscattered electrons
                           'p1EInf': 0.02,
                           'p1Ehat': 0.496,
                           'eEHat': 0.,
                           'w': 60.86,
                           'p': 1.,
                           'e1': 0.26,
                           'e2': 2.,
                           'sigmaE': 2.,
                           # Parameters for rediffused electrons
                           'p1RInf': 0.2,
                           'eR': 0.041,
                           'r': 0.104,
                           'q': 0.5,
                           'r1': 0.26,
                           'r2': 2.,
                           # Parameters for true secondaries
                           'deltaTSHat': 1.8848,
                           'eHat0': 332.,
                           's': 1.35,
                           't1': 0.5,  # t1 and t2 based on taylor expansion
                           't2': 1.,   # of PyECLOUD formula for E_max(theta)
                           't3': 0.7,
                           't4': 1.,
                           }
furman_pivi_surface = {'conserve_energy': False,
                       'exclude_rediffused': False,
                       'choice': 'poisson',
                       'M_cut': 10,
                       'p_n': np.array([2.5, 3.3, 2.5, 2.5, 2.8, 1.3, 1.5, 1.5, 1.5, 1.5]),
                       'eps_n': np.array([1.5, 1.75, 1., 3.75, 8.5, 11.5, 2.5, 3., 2.5, 3.]),
                       'p1EInf': 0.02,
                       'p1Ehat': 0.496,
                       'eEHat': 0.,
                       'w': 60.86,
                       'p': 1.,
                       'e1': 0.26,
                       'e2': 2.,
                       'sigmaE': 2.,
                       'p1RInf': 0.2,
                       'eR': 0.041,
                       'r': 0.104,
                       'q': 0.5,
                       'r1': 0.26,
                       'r2': 2.,
                       'deltaTSHat': 1.8848,
                       'eHat0': 276.8,
                       's': 1.54,
                       't1': 0.66,
                       't2': 0.8,
                       't3': 0.7,
                       't4': 1.,
                       }

sey_mod = fp.SEY_model_furman_pivi(E_th=35., sigmafit=1.0828, mufit=1.6636, secondary_angle_distribution='cosine_3D',
                                   switch_no_increase_energy=0, thresh_low_energy=-1,
                                   furman_pivi_surface=furman_pivi_surface_LHC)

chamb = ellip_cham_geom_object(1., 1., flag_verbose_file=False)
impact_management_object = impact_management(chamb=chamb, sey_mod=sey_mod, Dx_hist=.1, scrub_en_th=25.,
                                             Nbin_En_hist=100, En_hist_max=3000, flag_seg=False,
                                             cos_angle_width=0.05, flag_cos_angle_hist=True)

cos_theta_test = np.linspace(.1, 1., 10)
E_0_single = 35.
n_rep = int(1e4)
E_impact_eV_test = np.array([E_0_single] * n_rep)
alpha = 0.9

dists = impact_management_object.extract_energy_distributions(n_rep, E_impact_eV_test, cos_theta_test, mass=me)

plt.close('all')
ms.mystyle_arial()

fig1 = plt.figure(1, figsize=(3 * 8, 2 * 8))
fig1.set_facecolor('w')
sp1 = fig1.add_subplot(2, 2, 1)
sp2 = fig1.add_subplot(2, 2, 2)
sp3 = fig1.add_subplot(2, 2, 3)
sp4 = fig1.add_subplot(2, 2, 4)

for i_ct, ct in enumerate(cos_theta_test):
    thiscol = ms.colorprog(i_ct, len(cos_theta_test))
    label = 'costheta=%.2f' % ct
    sp1.hist(dists['true'][i_ct], bins=60, color=thiscol, label=label, alpha=alpha, density=True)
    sp2.hist(dists['elast'][i_ct], bins=30, color=thiscol, label=label, alpha=alpha, density=True)
    sp3.hist(dists['rediff'][i_ct], bins=30, color=thiscol, label=label, alpha=alpha, density=True)
    sp4.hist(dists['absorb'][i_ct], bins=30, color=thiscol, label=label, alpha=alpha, density=True)

linewid = 3
sp2.plot(0, 0, 'k', label='Model PDF', linewidth=linewid)
sp2.legend(loc='best', prop={'size': 14})
sz = 35
sp1.set_ylabel('True secondaries', fontsize=sz)
sp2.set_ylabel('Elastic', fontsize=sz)
sp3.set_ylabel('Rediffused', fontsize=sz)
sp4.set_ylabel('Absorbed', fontsize=sz)

# Compare with model
test_obj = sey_mod
E_0 = np.array([E_0_single] * int(1e5))
energy = np.linspace(0.001, E_0_single, num=int(1e5))
# Rediffused
prob_density_r = test_obj.rediffused_energy_PDF(energy=energy, E_0=E_0)
sp3.plot(energy, prob_density_r, 'k', label='PDF', linewidth=linewid)
# Backscaterred
prob_density_e = test_obj.backscattered_energy_PDF(energy, E_0)
sp2.plot(energy, prob_density_e, 'k', label='PDF', linewidth=linewid)
# True secondaries
delta_e, delta_r, delta_ts = test_obj.yield_fun_furman_pivi(E=E_0_single, costheta=1.)
delta_ts_prime = delta_ts / (1 - delta_e - delta_r)  # delta_ts^prime in FP paper, eq. (39)
prob_density_ts = test_obj.average_true_sec_energy_PDF(delta_ts=delta_ts_prime, E_0=E_0_single, energy=energy)
sp1.plot(energy, prob_density_ts, 'k', label='PDF of true secondary electrons (average)', linewidth=linewid)

for sp in [sp1, sp2, sp3, sp4]:
    sp.grid('on')
    sp.set_xlabel('Electron energy [eV]')

plt.subplots_adjust(right=0.99, left=.075)

plt.suptitle('Energy distribution extraction tests: Furman-Pivi model', fontsize=30)

plt.figure(2, figsize=(12, 9), facecolor='white')
for M in np.arange(1, sey_mod.M_cut + 1, 1):
    prob_density_ts = test_obj.true_sec_energy_PDF(delta_ts=delta_ts_prime, nn=M, E_0=E_0_single, energy=energy)
    plt.plot(energy, prob_density_ts[0], label='n: %i' % M, linewidth=linewid)
plt.legend()
plt.title('Energy distribution PDFs for secondary electron energies \n in the Furman-Pivi model', fontsize=sz - 10)
plt.xlabel('Energy [eV]', fontsize=sz - 10)
plt.ylabel('Normalised energy spectrum', fontsize=sz - 10)
plt.grid(.3)

plt.show()

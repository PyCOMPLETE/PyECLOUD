import sys
if '../..' not in sys.path:
    sys.path.append('../..')
import numpy as np
import matplotlib.pyplot as plt
import mystyle as ms
from scipy.optimize import curve_fit
import sec_emission_model_furman_pivi as fp
from impact_management_class import impact_management
from geom_impact_ellip import ellip_cham_geom_object

me = 9.10938356e-31
alpha = 0.9

furman_pivi_surface_LHC = {
    'use_modified_sigmaE': False,
    'use_ECLOUD_theta0_dependence': False,
    'use_ECLOUD_energy': False,
    'conserve_energy': False,
    'exclude_rediffused': True,
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
    't4': 1.}

flag_costheta_delta_scale = True
flag_costheta_Emax_shift = True

sey_mod_FP = fp.SEY_model_furman_pivi(furman_pivi_surface_LHC,
                                      E_th=None, sigmafit=None, mufit=None,
                                      switch_no_increase_energy=0, thresh_low_energy=None,
                                      secondary_angle_distribution='cosine_3D',
                                      flag_costheta_Emax_shift=flag_costheta_Emax_shift,
                                      flag_costheta_delta_scale=flag_costheta_delta_scale)

chamb = ellip_cham_geom_object(1., 1., flag_verbose_file=False)
impact_management_object = impact_management(chamb=chamb, sey_mod=sey_mod_FP, Dx_hist=.1, scrub_en_th=25.,
                                             Nbin_En_hist=100, En_hist_max=3000, flag_seg=False,
                                             cos_angle_width=0.05, flag_cos_angle_hist=True)

plt.close('all')
ms.mystyle(20)
linewid = 3
fontsz = 25
alph = 0.5


# Elastic
def del_elas_ECLOUD(energy, R_0=0.7, E_max=332., E_0=150.):
    del_elas = R_0 * ((np.sqrt(energy) - np.sqrt(energy + E_0)) / (np.sqrt(energy) + np.sqrt(energy + E_0)))**2
    return del_elas


def del_elas_FP(energy, p1EInf, p1Ehat, ww, pp, p0=[1., 1., 1., 1.]):
    eEHat = 0
    exp_factor = -(np.abs(energy - eEHat) / ww)**pp / pp
    delta_e0 = p1EInf + (p1Ehat - p1EInf) * np.exp(exp_factor)
    return delta_e0


energy = np.linspace(0., 800, num=int(1e3))
del_e_ECLOUD = del_elas_ECLOUD(energy)

popt, pcov = curve_fit(del_elas_FP, energy, del_e_ECLOUD)

fitted = del_elas_FP(energy, *popt)

fig = plt.figure(figsize=(24, 9), facecolor='white')
fig.suptitle('Fitting Furman-Pivi model to ECLOUD model', fontsize=fontsz)
sp1 = fig.add_subplot(1, 2, 1)
sp1.plot(energy, del_e_ECLOUD, color='k', label='ECLOUD', linewidth=linewid)
sp1.plot(energy, fitted, '--', color='r', label='Furman-Pivi', linewidth=linewid)
sp1.legend()
sp1.set_ylabel(r'$\delta_{elas}$', fontsize=fontsz)
sp1.set_xlabel('Impacting electron energy [eV]', fontsize=fontsz)
sp1.grid(alpha=alph)
sp1.text(310, 0.3, 'Fitted parameters: \np1EInf: %f \np1EHat: %f \nW: %f \np: %f' % (popt[0], popt[1], popt[2], popt[3]))


# True secondary
def del_true_ECLOUD(energy, cos_theta, del_max, s=1.35, E_max=332.):
    E_max_tilde = E_max * (1. + 0.7 * (1. - cos_theta))
    x = energy / E_max_tilde
    del_true = del_max * s * x / (s - 1 + x**s)
    angular_factor = np.exp((1. - cos_theta) / 2.)
    return del_true * angular_factor


def _D(x, s):
    return s * x / (s - 1 + x**s)


energy = np.linspace(0.00001, 3000, num=int(1e3))
del_max = 1.6
cos_theta = .7


def delta_ts(energy, cos_theta, t1, t2):
    s = 1.35
    eHat0 = 332.
    deltaTSHat = del_max
    t3 = 0.7
    t4 = 1.
    eHat = eHat0 * (1. + t3 * (1. - cos_theta**t4))
    delta_ts0 = deltaTSHat * _D(energy / eHat, s)
    angular_factor = 1. + t1 * (1. - cos_theta**t2)

    return delta_ts0 * angular_factor


def fit_func(energy, t1, t2):
    return delta_ts(energy, cos_theta, t1, t2)


del_ts_ECLOUD = del_true_ECLOUD(energy, cos_theta, del_max)

popt_true, pcov_true = curve_fit(fit_func, energy, del_ts_ECLOUD)

fitted_true = fit_func(energy, *popt_true)

taylor_true = delta_ts(energy, cos_theta, t1=0.675766, t2=0.767523)

sp2 = fig.add_subplot(1, 2, 2)
sp2.plot(energy, del_ts_ECLOUD, color='k', label='ECLOUD', linewidth=linewid)
sp2.plot(energy, fitted_true, '--', color='r', label='Furman-Pivi', linewidth=linewid)
sp2.plot(energy, taylor_true, linestyle='-.', color='c', label='Furman-Pivi, Taylor', linewidth=linewid)
sp2.legend()
sp2.set_ylabel(r'$\delta_{true}$', fontsize=fontsz)
sp2.set_xlabel('Impacting electron energy [eV]', fontsize=fontsz)
sp2.grid(alpha=alph)
sp2.text(310, 0.3, 'Fitted parameters: \nt1: %f \nt2: %f' % (popt_true[0], popt_true[1]))
fig.subplots_adjust(left=0.06, right=0.95)

print("NOTE: The angular dependence of del_true is not the same, t1 and t2 can't be chosen such that there is good agreement for all impacting angles.")

plt.show()

import sys
if '../..' not in sys.path:
    sys.path.append('../..')
import numpy as np
import matplotlib.pyplot as plt
import mystyle as ms
from scipy.optimize import curve_fit
from scipy.special import binom
from scipy.misc import factorial
from scipy.integrate import simps
import sec_emission_model_furman_pivi as fp


plt.close('all')
ms.mystyle(20)
linewid = 2
fontsz = 20
alph = 0.5


# Elastic component
def del_elas_ECLOUD(energy, R_0=0.7, E_max=332., E_0=150.):
    del_elas = R_0 * ((np.sqrt(energy) - np.sqrt(energy + E_0)) / (np.sqrt(energy) + np.sqrt(energy + E_0)))**2
    return del_elas


def del_elas_FP(energy, p1EInf, p1Ehat, ww, pp):
    eEHat = 0
    exp_factor = -(np.abs(energy - eEHat) / ww)**pp / pp
    delta_e0 = p1EInf + (p1Ehat - p1EInf) * np.exp(exp_factor)
    return delta_e0


energy = np.linspace(0., 800, num=int(1e3))
del_e_ECLOUD = del_elas_ECLOUD(energy)

popt, pcov = curve_fit(del_elas_FP, energy, del_e_ECLOUD)

fitted = del_elas_FP(energy, *popt)

fig = plt.figure(figsize=(20, 18), facecolor='white')
fig.suptitle('Fitting Furman-Pivi model to ECLOUD model')
sp1 = fig.add_subplot(2, 2, 1)
sp1.plot(energy, del_e_ECLOUD, color='k', label='ECLOUD', linewidth=linewid)
sp1.plot(energy, fitted, '--', color='r', label='Furman-Pivi', linewidth=linewid)
sp1.legend()
sp1.set_ylabel(r'$\delta_{elas}$', fontsize=fontsz)
sp1.set_xlabel('Impacting electron energy [eV]', fontsize=fontsz)
sp1.grid(alpha=alph)
sp1.text(310, 0.3, 'Fitted parameters: \np1EInf: %f \np1EHat: %f \nW: %f \np: %f' % (popt[0], popt[1], popt[2], popt[3]))


# True secondary component
def del_true_ECLOUD(energy, cos_theta, del_max, s=1.35, E_max=332.):
    E_max_tilde = E_max * (1. + 0.7 * (1. - cos_theta))
    x = energy / E_max_tilde
    del_true = del_max * s * x / (s - 1 + x**s)
    angular_factor = np.exp((1. - cos_theta) / 2.)
    return del_true * angular_factor


def _D(x, s):
    return s * x / (s - 1 + x**s)


energy = np.linspace(0., 3000, num=int(1e3))
del_max = 1.8848
cos_theta = 0.9


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

sp2 = fig.add_subplot(2, 2, 2)
sp2.plot(energy, del_ts_ECLOUD, color='k', label='ECLOUD', linewidth=linewid)
sp2.plot(energy, fitted_true, '--', color='r', label='Furman-Pivi', linewidth=linewid)
sp2.legend()
sp2.set_ylabel(r'$\delta_{true}$', fontsize=fontsz)
sp2.set_xlabel('Impacting electron energy [eV]', fontsize=fontsz)
sp2.grid(alpha=alph)
sp2.text(310, 0.3, 'Fitted parameters: \nt1: %f \nt2: %f' % (popt_true[0], popt_true[1]))
fig.subplots_adjust(left=0.06, right=0.95)


# Energy distributions
M_cut = 10
choice = 'poisson'


def true_sec_energy_PDF(delta_ts, nn, E_0, energy):
    """Gives the PDF for true secondary electrons."""
    if nn == 0:
        raise ValueError('nn = 0, you cannot emit zero electrons.')
    p_n = np.array([2.5, 3.3, 2.5, 2.5, 2.8, 1.3, 1.5, 1.5, 1.5, 1.5])
    eps_n = np.array([1.5, 1.75, 1., 3.75, 8.5, 11.5, 2.5, 3., 2.5, 3.])

    nn_all = np.arange(0, M_cut + 1, 1)

    if choice == 'poisson':
        P_n_ts = np.squeeze(delta_ts**nn_all / factorial(nn_all) * np.exp(-delta_ts))
    elif choice == 'binomial':
        p = delta_ts / M_cut
        P_n_ts = np.squeeze(binom(M_cut, nn) * (p)**nn_all * (1 - p)**(M_cut - nn_all))
    else:
        raise ValueError('choice must be either \'poisson\' or \'binomial\'')

    P_n_ts = P_n_ts / np.sum(P_n_ts)
    P_n_ts_return = P_n_ts[int(nn)]
    eps_curr = eps_n[int(nn - 1)]
    p_n_curr = p_n[int(nn - 1)]

    if E_0 == 0:
        F_n = 0
    else:
        F_n = 1
    f_n_ts = F_n * energy**(p_n_curr - 1) * np.exp(-energy / eps_curr)
    area = simps(f_n_ts, energy)
    f_n_ts = f_n_ts / area  # normalisation

    return f_n_ts, P_n_ts_return


def average_true_sec_energy_PDF(delta_ts, E_0, energy):
    nns = np.arange(1, self.M_cut + 1, 1)
    average_f_n_ts = np.zeros_like(energy)
    for ii in nns:
        f_n_ts, P_n_ts = self.true_sec_energy_PDF(delta_ts=delta_ts, nn=ii, E_0=E_0, energy=energy)
        average_f_n_ts = average_f_n_ts + f_n_ts * P_n_ts * ii
    area = simps(average_f_n_ts, energy)
    return average_f_n_ts / area


def hilleret_energy(energy):
    sigmafit = 1.0828
    mufit = 1.6636
    return 1. / (energy * sigmafit * np.sqrt(2 * np.pi)) * np.exp(-(np.log(energy) - mufit)**2 / (2 * sigmafit**2))


def simple_pdf(energy, pp, eps):
    pdf = energy**(pp - 1) * np.exp(-energy / eps)
    area = simps(pdf, energy)
    return pdf / area


ene_hilleret = hilleret_energy(energy)

popt_ene_true, pcov_ene_true = curve_fit(simple_pdf, energy, ene_hilleret)
fitted_ene_true = simple_pdf(energy, *popt_ene_true)

sp3 = fig.add_subplot(2, 2, 3)
sp3.plot(energy, ene_hilleret, color='k', label='ECLOUD', linewidth=linewid)
sp3.plot(energy, fitted_ene_true, '--', color='r', label='Furman-Pivi', linewidth=linewid)
sp3.legend()
sp3.set_ylabel('Energy distribution', fontsize=fontsz)
sp3.set_xlabel('Emitted electron energy [eV]', fontsize=fontsz)
sp3.grid(alpha=alph)
sp3.text(310, 0.3, 'Fitted parameters: \np: %f \neps: %f' % (popt_ene_true[0], popt_ene_true[1]))



# furman_pivi_surface_tweaked = {'exclude_rediffused': False,
#                                'choice': 'poisson',
#                                'M_cut': 10,
#                                'p_n': np.array([2.5, 3.3, 2.5, 2.5, 2.8, 1.3, 1.5, 1.5, 1.5, 1.5]),
#                                'eps_n': np.array([1.5, 1.75, 1., 3.75, 8.5, 11.5, 2.5, 3., 2.5, 3.]),
#                                # Parameters for backscattered electrons
#                                'p1EInf': 0.002158,  # Changed this
#                                'p1Ehat': 0.709633,  # Changed this
#                                'eEHat': 0.,
#                                'w': 46.028959,  # Changed this
#                                'p': 0.468907,  # Changed this
#                                'e1': 0.,  # Changed this
#                                'e2': 2.,
#                                'sigmaE': 2.,
#                                # Parameters for rediffused electrons
#                                'p1RInf': 0.2,
#                                'eR': 0.041,
#                                'r': 0.104,
#                                'q': 0.5,
#                                'r1': 0.26,
#                                'r2': 2.,
#                                # Parameters for true secondaries
#                                'deltaTSHat': 1.8848,
#                                'eHat0': 332.,
#                                's': 1.35,
#                                't1': 0.727814,  # Changed this
#                                't2': 0.69333,  # Changed this
#                                't3': 0.7,
#                                't4': 1.,
#                                }
#
# sey_mod = fp.SEY_model_furman_pivi(E_th=35., sigmafit=1.0828, mufit=1.6636, secondary_angle_distribution='cosine_3D',
#                                    switch_no_increase_energy=0, thresh_low_energy=-1,
#                                    furman_pivi_surface=furman_pivi_surface_tweaked)


plt.show()

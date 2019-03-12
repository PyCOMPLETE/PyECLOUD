import sys
if '../..' not in sys.path:
    sys.path.append('../..')
import numpy as np
import matplotlib.pyplot as plt
import mystyle as ms
from scipy.optimize import curve_fit
from scipy.misc import factorial
from scipy.integrate import simps
import scipy


plt.close('all')
ms.mystyle(20)
linewid = 2
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


def pdf_elas_FP(energy, E_imp, sigma_e=2.):
    return 2 * np.exp(-(energy - E_imp)**2 / (2 * sigma_e**2)) / (np.sqrt(2 * np.pi) * sigma_e * scipy.special.erf(E_imp / np.sqrt(2) / sigma_e))


energy = np.linspace(0., 800, num=int(1e3))
del_e_ECLOUD = del_elas_ECLOUD(energy)

popt, pcov = curve_fit(del_elas_FP, energy, del_e_ECLOUD)

fitted = del_elas_FP(energy, *popt)

fig = plt.figure(figsize=(20, 18), facecolor='white')
fig.suptitle('Fitting Furman-Pivi model to ECLOUD model', fontsize=fontsz)
sp1 = fig.add_subplot(3, 2, 1)
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

sp2 = fig.add_subplot(3, 2, 2)
sp2.plot(energy, del_ts_ECLOUD, color='k', label='ECLOUD', linewidth=linewid)
sp2.plot(energy, fitted_true, '--', color='r', label='Furman-Pivi', linewidth=linewid)
sp2.legend()
sp2.set_ylabel(r'$\delta_{true}$', fontsize=fontsz)
sp2.set_xlabel('Impacting electron energy [eV]', fontsize=fontsz)
sp2.grid(alpha=alph)
sp2.text(310, 0.3, 'Fitted parameters: \nt1: %f \nt2: %f' % (popt_true[0], popt_true[1]))
fig.subplots_adjust(left=0.06, right=0.95)


def hilleret_energy(energy):
    sigmafit = 1.0828
    mufit = 1.6636
    return 1. / (energy * sigmafit * np.sqrt(2 * np.pi)) * np.exp(-(np.log(energy) - mufit)**2 / (2 * sigmafit**2))


def simple_pdf(energy, pp, eps):
    pdf = energy**(pp - 1) * np.exp(-energy / eps)
    area = simps(pdf, energy)
    return pdf / area


energy = np.linspace(0.0000001, 35., num=int(1e3))
ene_hilleret = hilleret_energy(energy)

popt_ene_true, pcov_ene_true = curve_fit(simple_pdf, energy, ene_hilleret, p0=[1., 1.])
fitted_ene_true = simple_pdf(energy, *popt_ene_true)

sp3 = fig.add_subplot(3, 2, 3)
sp3.plot(energy, ene_hilleret, color='k', label='ECLOUD', linewidth=linewid)
sp3.plot(energy, fitted_ene_true, '--', color='r', label='Furman-Pivi', linewidth=linewid)
sp3.legend()
sp3.set_ylabel('Energy distribution', fontsize=fontsz)
sp3.set_xlabel('Emitted electron energy [eV]', fontsize=fontsz)
sp3.grid(alpha=alph)
sp3.text(20, 0.04, 'Fitted parameters: \np: %f \neps: %f' % (popt_ene_true[0], popt_ene_true[1]))


M_cut = 10


def true_sec_energy_PDF(nn, energy, p_n, eps_n):
    """The PDF for true secondary electrons."""
    if nn == 0:
        raise ValueError('nn = 0, you cannot emit zero electrons.')

    delta_ts = 1.88

    nn_all = np.arange(0, M_cut + 1, 1)

    P_n_ts = np.squeeze(delta_ts**nn_all / factorial(nn_all) * np.exp(-delta_ts))

    P_n_ts = P_n_ts / np.sum(P_n_ts)
    P_n_ts_return = P_n_ts[int(nn)]
    eps_curr = eps_n[int(nn - 1)]
    p_n_curr = p_n[int(nn - 1)]

    F_n = 1
    f_n_ts = F_n * energy**(p_n_curr - 1) * np.exp(-energy / eps_curr)
    area = scipy.integrate.simps(f_n_ts, energy)
    f_n_ts = f_n_ts / area  # normalisation

    return f_n_ts, P_n_ts_return


def average_true_sec_energy_PDF(energy, p_1, p2, p3, p4, p5, p6, p7, p8, p9, p10,
                                eps_1, eps2, eps3, eps4, eps5, eps6, eps7, eps8, eps9, eps10):

    p_n = np.array([p_1, p2, p3, p4, p5, p6, p7, p8, p9, p10])
    eps_n = np.array([eps_1, eps2, eps3, eps4, eps5, eps6, eps7, eps8, eps9, eps10])

    nns = np.arange(1, M_cut + 1, 1)
    average_f_n_ts = np.zeros_like(energy)
    for ii in nns:
        f_n_ts, P_n_ts = true_sec_energy_PDF(nn=ii, energy=energy, p_n=p_n, eps_n=eps_n)
        average_f_n_ts = average_f_n_ts + f_n_ts * P_n_ts * ii
    area = scipy.integrate.simps(average_f_n_ts, energy)
    return average_f_n_ts / area


energy = np.linspace(0.00001, 35., num=int(1e3))
ene_hilleret = hilleret_energy(energy)

popt_ene_true, pcov_ene_true = curve_fit(average_true_sec_energy_PDF, energy, ene_hilleret,
                                         p0=[1.] * 20, bounds=(0.01, 100.))
fitted_ene_true = average_true_sec_energy_PDF(energy, *popt_ene_true)

sp4 = fig.add_subplot(3, 2, 4)
sp4.plot(energy, ene_hilleret, color='k', label='ECLOUD', linewidth=linewid)
sp4.plot(energy, fitted_ene_true, '--', color='r', label='Furman-Pivi', linewidth=linewid)
sp4.legend()
sp4.set_ylabel('Energy distribution', fontsize=fontsz)
sp4.set_xlabel('Emitted electron energy [eV]', fontsize=fontsz)
sp4.grid(alpha=alph)
sp4.text(4, 0.068, 'Fitted parameters: \np_n: %s \n %s \n %s \nep_n: %s \n %s \n %s' % (str(popt_ene_true[0:4]),
                                                                                        str(popt_ene_true[4:7]),
                                                                                        str(popt_ene_true[7:10]),
                                                                                        str(popt_ene_true[10:14]),
                                                                                        str(popt_ene_true[14:17]),
                                                                                        str(popt_ene_true[17:])))
print('Fitted parameters for average true energy PDF:')
print('p_n:')
print(popt_ene_true[0:10])
print('eps_n:')
print(popt_ene_true[10:])


# Comparison with other impactin energies
# energy = np.linspace(0.00001, 100., num=int(1e3))
# ene_hilleret = hilleret_energy(energy)
#
# fitted_ene_true = average_true_sec_energy_PDF(energy, *popt_ene_true)
#
# sp5 = fig.add_subplot(3, 2, 5)
# sp5.plot(energy, ene_hilleret, color='k', label='ECLOUD', linewidth=linewid)
# sp5.plot(energy, fitted_ene_true, '--', color='r', label='Furman-Pivi', linewidth=linewid)
# sp5.legend()
# sp5.set_ylabel('Energy distribution', fontsize=fontsz)
# sp5.set_xlabel('Emitted electron energy [eV]', fontsize=fontsz)
# sp5.grid(alpha=alph)
#
# energy = np.linspace(0.00001, 10., num=int(1e3))
# ene_hilleret = hilleret_energy(energy)
#
# fitted_ene_true = average_true_sec_energy_PDF(energy, *popt_ene_true)
#
# sp6 = fig.add_subplot(3, 2, 6)
# sp6.plot(energy, ene_hilleret, color='k', label='ECLOUD', linewidth=linewid)
# sp6.plot(energy, fitted_ene_true, '--', color='r', label='Furman-Pivi', linewidth=linewid)
# sp6.legend()
# sp6.set_ylabel('Energy distribution', fontsize=fontsz)
# sp6.set_xlabel('Emitted electron energy [eV]', fontsize=fontsz)
# sp6.grid(alpha=alph)


# energy = np.linspace(0.0000001, 35., num=int(1e3))
# ene_FP = pdf_elas_FP(energy, E_imp=35., sigma_e=.01)
#
# popt_ene_elas, pcov_ene_elas = curve_fit(pdf_elas_FP, energy, ene_FP, p0=[1., 1.])
# fitted_ene_elas = pdf_elas_FP(energy, *popt_ene_elas)
#
# sp4 = fig.add_subplot(2, 2, 4)
# sp4.plot(energy, ene_FP, color='k', label='ECLOUD', linewidth=linewid)
# sp4.plot(energy, fitted_ene_elas, '--', color='r', label='Furman-Pivi', linewidth=linewid)
# sp4.legend()
# sp4.set_ylabel('Energy distribution', fontsize=fontsz)
# sp4.set_xlabel('Emitted electron energy [eV]', fontsize=fontsz)
# sp4.grid(alpha=alph)
# sp4.text(20, 0.04, 'Fitted parameters: \nE_0: %f \n$\sigma_e$: %f' % (popt_ene_elas[0], popt_ene_elas[1]))

plt.show()

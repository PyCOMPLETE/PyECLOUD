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
import sec_emission_model_furman_pivi as fp
from impact_management_class import impact_management
from geom_impact_ellip import ellip_cham_geom_object

me = 9.10938356e-31
alpha = 0.9

furman_pivi_surface_LHC = {'exclude_rediffused': True,
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


sey_mod_FP = fp.SEY_model_furman_pivi(furman_pivi_surface_LHC,
                                      E_th=None, sigmafit=None, mufit=None,
                                      switch_no_increase_energy=0, thresh_low_energy=None,
                                      secondary_angle_distribution='cosine_3D')

chamb = ellip_cham_geom_object(1., 1., flag_verbose_file=False)
impact_management_object = impact_management(chamb=chamb, sey_mod=sey_mod_FP, Dx_hist=.1, scrub_en_th=25.,
                                             Nbin_En_hist=100, En_hist_max=3000, flag_seg=False,
                                             cos_angle_width=0.05, flag_cos_angle_hist=True)

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
    hilleret = 1. / (energy * sigmafit * np.sqrt(2 * np.pi)) * np.exp(-(np.log(energy) - mufit)**2 / (2 * sigmafit**2))
    area = simps(hilleret, energy)
    return hilleret / area

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


E_imp_used_for_fit = 100.
delta_e, _, delta_ts = sey_mod_FP.yield_fun_furman_pivi(E=E_imp_used_for_fit, costheta=1.)
delta_ts_prime = delta_ts / (1 - delta_e)


def true_sec_energy_PDF(nn, energy, p_n, eps_n):
    """The PDF for true secondary electrons."""
    if nn == 0:
        raise ValueError('nn = 0, you cannot emit zero electrons.')

    delta_ts = delta_ts_prime

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


energy = np.linspace(0.00001, E_imp_used_for_fit, num=int(1e3))
ene_hilleret = hilleret_energy(energy)

bounds = ([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
           1., 1., 1., 1., 1., 1., 1., 1., 1., 1.],
          [4., 4., 4., 4., 30., 30., 30., 30., 30., 30.,
           4., 4., 4., 5., 30., 30., 30., 30., 30., 30.])

popt_ene_true, pcov_ene_true = curve_fit(average_true_sec_energy_PDF, energy, ene_hilleret,
                                         bounds=(1., 30.))
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
energy = np.linspace(0.00001, 100., num=int(1e3))
ene_hilleret = hilleret_energy(energy)

delta_e, _, delta_ts = sey_mod_FP.yield_fun_furman_pivi(E=energy[-1], costheta=1.)
delta_ts_prime = delta_ts / (1 - delta_e)


def true_sec_energy_PDF(nn, energy, p_n, eps_n):
    """The PDF for true secondary electrons."""
    if nn == 0:
        raise ValueError('nn = 0, you cannot emit zero electrons.')

    delta_ts = delta_ts_prime

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

fitted_ene_true = average_true_sec_energy_PDF(energy, *popt_ene_true)


sp5 = fig.add_subplot(3, 2, 5)
sp5.plot(energy, ene_hilleret, color='k', label='ECLOUD', linewidth=linewid)
sp5.plot(energy, fitted_ene_true, '--', color='r', label='Furman-Pivi', linewidth=linewid)
sp5.legend()
sp5.set_ylabel('Energy distribution', fontsize=fontsz)
sp5.set_xlabel('Emitted electron energy [eV]', fontsize=fontsz)
sp5.grid(alpha=alph)

energy = np.linspace(0.00001, 20., num=int(1e3))
ene_hilleret = hilleret_energy(energy)

delta_e, _, delta_ts = sey_mod_FP.yield_fun_furman_pivi(E=energy[-1], costheta=1.)
delta_ts_prime = delta_ts / (1 - delta_e)


def true_sec_energy_PDF(nn, energy, p_n, eps_n):
    """The PDF for true secondary electrons."""
    if nn == 0:
        raise ValueError('nn = 0, you cannot emit zero electrons.')

    delta_ts = delta_ts_prime

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


fitted_ene_true = average_true_sec_energy_PDF(energy, *popt_ene_true)

sp6 = fig.add_subplot(3, 2, 6)
sp6.plot(energy, ene_hilleret, color='k', label='ECLOUD', linewidth=linewid)
sp6.plot(energy, fitted_ene_true, '--', color='r', label='Furman-Pivi', linewidth=linewid)
sp6.legend()
sp6.set_ylabel('Energy distribution', fontsize=fontsz)
sp6.set_xlabel('Emitted electron energy [eV]', fontsize=fontsz)
sp6.grid(alpha=alph)


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

plt.figure(2)
E_0s = np.array([5., 10., 20., 35., 50., 75., 100.])
for i_E, E_0_curr in enumerate(E_0s):
    thiscol = ms.colorprog(i_E, len(E_0s))
    energy = np.linspace(0.00001, E_0_curr, num=int(1e3))
    ene_hilleret = hilleret_energy(energy)
    delta_e, _, delta_ts = sey_mod_FP.yield_fun_furman_pivi(E=E_0_curr, costheta=1.)
    delta_ts_prime = delta_ts / (1 - delta_e)

    dists = impact_management_object.extract_energy_distributions(n_rep=int(1e5), E_impact_eV_test=np.array([E_0_curr] * int(1e5)), cos_theta_test=[1.], mass=me)


    # def true_sec_energy_PDF(nn, energy, p_n, eps_n):
    #     """The PDF for true secondary electrons."""
    #     if nn == 0:
    #         raise ValueError('nn = 0, you cannot emit zero electrons.')
    #
    #     delta_ts = delta_ts_prime
    #
    #     nn_all = np.arange(0, M_cut + 1, 1)
    #
    #     P_n_ts = np.squeeze(delta_ts**nn_all / factorial(nn_all) * np.exp(-delta_ts))
    #
    #     P_n_ts = P_n_ts / np.sum(P_n_ts)
    #     P_n_ts_return = P_n_ts[int(nn)]
    #     eps_curr = eps_n[int(nn - 1)]
    #     p_n_curr = p_n[int(nn - 1)]
    #
    #     F_n = 1
    #     f_n_ts = F_n * energy**(p_n_curr - 1) * np.exp(-energy / eps_curr)
    #     area = scipy.integrate.simps(f_n_ts, energy)
    #     f_n_ts = f_n_ts / area  # normalisation
    #
    #     return f_n_ts, P_n_ts_return
    #
    #
    # def average_true_sec_energy_PDF(energy, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,
    #                                 eps1, eps2, eps3, eps4, eps5, eps6, eps7, eps8, eps9, eps10):
    #
    #     p_n = np.array([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10])
    #     eps_n = np.array([eps1, eps2, eps3, eps4, eps5, eps6, eps7, eps8, eps9, eps10])
    #
    #     nns = np.arange(1, M_cut + 1, 1)
    #     average_f_n_ts = np.zeros_like(energy)
    #     for ii in nns:
    #         f_n_ts, P_n_ts = true_sec_energy_PDF(nn=ii, energy=energy, p_n=p_n, eps_n=eps_n)
    #         average_f_n_ts = average_f_n_ts + f_n_ts * P_n_ts * ii
    #     area = scipy.integrate.simps(average_f_n_ts, energy)
    #     return average_f_n_ts / area
    #
    # p1 = 2.5
    # p2 = 3.3
    # p3 = 2.5
    # p4 = 2.5
    # p5 = 2.8
    # p6 = 1.3
    # p7 = 1.5
    # p8 = 1.5
    # p9 = 1.5
    # p10 = 1.5
    # eps1 = 1.5
    # eps2 = 1.75
    # eps3 = 1.
    # eps4 = 3.75
    # eps5 = 8.5
    # eps6 = 11.5
    # eps7 = 2.5
    # eps8 = 3.0
    # eps9 = 2.5
    # eps10 = 3.0
    #
    # pdf = average_true_sec_energy_PDF(energy, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,
    #                                   eps1, eps2, eps3, eps4, eps5, eps6, eps7, eps8, eps9, eps10)
    # plt.plot(energy, pdf, color=thiscol, linewidth=linewid, label='$E_0 =$ %.0f' % E_0_curr)
    # pdf = sey_mod_FP.average_true_sec_energy_PDF(delta_ts=delta_ts_prime, E_0=E_imp_used_for_fit, energy=energy)

    pdf = average_true_sec_energy_PDF(energy, *popt_ene_true)
    plt.plot(energy, pdf, color=thiscol, linestyle='--', linewidth=linewid, label='$E_0 =$ %.0f' % E_0_curr)
    plt.plot(energy, ene_hilleret, color=thiscol, linestyle='-', linewidth=linewid, label='$E_0 =$ %.0f' % E_0_curr)
    plt.hist(dists['true'][0], color=thiscol, alpha=alpha, density=True)
plt.legend()


plt.show()

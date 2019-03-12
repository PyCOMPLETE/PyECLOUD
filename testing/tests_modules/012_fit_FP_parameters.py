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
M_cut = 10

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

sey_mod_initial = fp.SEY_model_furman_pivi(furman_pivi_surface_LHC,
                                           E_th=None, sigmafit=None, mufit=None,
                                           switch_no_increase_energy=0, thresh_low_energy=None,
                                           secondary_angle_distribution='cosine_3D')


def make_it_easier(*args):
    furman_pivi_surface_LHC = {'exclude_rediffused': True,
                               'choice': 'poisson',
                               'M_cut': M_cut,
                               'p_n': np.array(args[0:10]),
                               'eps_n': np.array(args[10:]),
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
    return sey_mod_FP, impact_management_object, furman_pivi_surface_LHC


plt.close('all')
ms.mystyle(20)
linewid = 2
fontsz = 25
alph = 0.5
E_imp_used_for_fit = 35.


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


delta_e, _, delta_ts = sey_mod_initial.yield_fun_furman_pivi(E=E_imp_used_for_fit, costheta=1.)
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


bounds = ([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
           1., 1., 1., 1., 1., 1., 1., 1., 1., 1.],
          [4., 4., 4., 4., 30., 30., 30., 30., 30., 30.,
           4., 4., 4., 5., 30., 30., 30., 30., 30., 30.])

energy = np.linspace(0.0001, E_imp_used_for_fit, 1000)
ene_hilleret = hilleret_energy(energy)
popt, pcov = curve_fit(average_true_sec_energy_PDF, energy, ene_hilleret, bounds=(1., 30.))


plt.figure(2)
E_0s = np.array([5., 10., 20., 35., 50., 75., 100.])
E_0s = np.array([20., 35., 50., 100.])

for i_E, E_0_curr in enumerate(E_0s):
    thiscol = ms.colorprog(i_E, len(E_0s))
    energy = np.linspace(0.00001, E_0_curr, num=int(1e3))
    ene_hilleret = hilleret_energy(energy)

    sey_mod_FP, impact_management_object, _ = make_it_easier(*popt)

    delta_e, _, delta_ts = sey_mod_FP.yield_fun_furman_pivi(E=E_0_curr, costheta=1.)
    delta_ts_prime = delta_ts / (1 - delta_e)

    dists = impact_management_object.extract_energy_distributions(n_rep=int(1e5), E_impact_eV_test=np.array([E_0_curr] * int(1e5)), cos_theta_test=[1], mass=me)

    pdf_from_module = sey_mod_FP.average_true_sec_energy_PDF(delta_ts=delta_ts_prime, E_0=E_0_curr, energy=energy)
    plt.plot(energy, pdf_from_module, 'k', linewidth=linewid)
    pdf = average_true_sec_energy_PDF(energy, *popt)
    plt.plot(energy, pdf, color=thiscol, linestyle='--', linewidth=linewid, label='FP, $E_0 =$ %.0f' % E_0_curr)
    plt.plot(energy, ene_hilleret, color=thiscol, linestyle='-', linewidth=linewid, label='ECLOUD, $E_0 =$ %.0f' % E_0_curr)
    # plt.hist(dists['true'][0], color=thiscol, alpha=alpha, density=True)
plt.legend()


plt.show()

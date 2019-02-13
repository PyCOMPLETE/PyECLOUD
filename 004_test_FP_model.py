import numpy as np
import matplotlib.pyplot as plt
import sec_emission_model_furman_pivi as fp
import mystyle as ms
import scipy

plt.close('all')
ms.mystyle(12)
linewid = 2

furman_pivi_surface_LHC = {'M': 10,
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
                           'eHat0': 322.,
                           's': 1.35,
                           't1': 0.66,
                           't2': 0.8,
                           't3': 0.7,
                           't4': 1.,
                           }
furman_pivi_surface = {'M': 10,
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
                       'eHat0': 322.,
                       's': 1.35,
                       't1': 0.66,
                       't2': 0.8,
                       't3': 0.7,
                       't4': 1.,
                       }

test_obj = fp.SEY_model_furman_pivi(E_th=35., sigmafit=1.0828, mufit=1.6636, secondary_angle_distribution='cosine_3D',
                                    switch_no_increase_energy=0, thresh_low_energy=-1,
                                    furman_pivi_surface=furman_pivi_surface_LHC)

qq = 0.5  # From FP paper
sigma_e = 2.
E_0_single = 100.
E_0 = np.array([E_0_single] * int(1e5))
energy = np.linspace(0.001, E_0_single, num=int(1e5))

alpha = 0.3
round_to_digits = 4

plt.close('all')
# Subplots
fig, ((ax1, ax2, ax4), (ax5, ax6, ax8)) = plt.subplots(2, 3, sharex=True, figsize=(1.8 * 12, 12), facecolor='w')
plt.suptitle('Furman-Pivi tests: Energy distributions', fontsize=30)

axlist = [ax1, ax2, ax4, ax5, ax6, ax8]
for ax in axlist:
    ax.grid(alpha=alpha)
    ax.set_xlabel('Energy [eV]')

# Backscattered
prob_density_e = test_obj.backscattered_energy_PDF(energy, E_0)
ax1.plot(energy, prob_density_e, label='PDF', linewidth=linewid)
ax5.plot(energy, test_obj.backscattered_energy_CDF(energy, E_0), label='CDF', linewidth=linewid)
ax1.set_title('Backscattered energy distribution')
area = scipy.integrate.simps(prob_density_e, energy)
area = round(area, round_to_digits)
ax1.text(E_0_single / 2., ax1.get_ylim()[1] / 2, 'Area = ' + str(area), fontsize=18)
ax1.legend()
ax5.legend()
ax1.hist(test_obj.get_energy_backscattered(E_0), density=True, bins=20)

# Rediffused
prob_density_r = test_obj.rediffused_energy_PDF(energy=energy, E_0=E_0)
ax2.plot(energy, prob_density_r, label='PDF', linewidth=linewid)
ax6.plot(energy, test_obj.rediffused_energy_CDF(energy=energy, E_0=E_0), label='CDF', linewidth=linewid)
ax2.set_title('Rediffused energy distribution')
area = scipy.integrate.simps(prob_density_r, energy)
area = round(area, round_to_digits)
ax2.text(E_0_single / 2., ax2.get_ylim()[1] / 2, 'Area = ' + str(area), fontsize=18)
ax2.hist(test_obj.get_energy_rediffused(E_0), density=True, bins=20)
ax2.legend()
ax6.legend()

# True secondary

delta_e, delta_r, delta_ts = test_obj._yield_fun_furman_pivi(E=E_0_single, costheta=1.)
delta_ts_prime = delta_ts / (1 - delta_e - delta_r)  # delta_ts^prime in FP paper, eq. (39)

nn = 2
prob_density_ts, pnts = test_obj.true_sec_energy_PDF(delta_ts=delta_ts, nn=nn, E_0=E_0_single, energy=energy)
ax4.plot(energy, prob_density_ts, label='PDF of true secondary electrons', linewidth=linewid)
ax4.set_title(r'True secondary energy distribution, $f_{%i,ts}$' % nn)
area = scipy.integrate.simps(prob_density_ts, energy)
area = round(area, round_to_digits)
ax4.text(E_0_single / 2., ax4.get_ylim()[1] / 2, 'Area = ' + str(area) + ', \nP_n_ts = ' + str(round(pnts, round_to_digits)), fontsize=18)
CDF, _ = test_obj.true_sec_energy_CDF(delta_ts=delta_ts, nn=nn, E_0=E_0_single, energy=energy)
ax8.plot(energy, CDF, label='CDF', linewidth=linewid)
ax4.hist(test_obj.get_energy_true_sec(delta_ts=delta_ts, nn=np.repeat(nn, len(E_0)), E_0=E_0), density=True, bins=30)
ax8.legend()

fig.subplots_adjust(left=0.05, right=0.95)

# SEY components
plt.figure(2, figsize=(16, 12), facecolor='w')
sp1 = plt.subplot(2, 2, 1)
sp2 = plt.subplot(2, 2, 2, sharex=sp1, sharey=sp1)
sp3 = plt.subplot(2, 2, 3, sharex=sp1, sharey=sp1)
sp4 = plt.subplot(2, 2, 4, sharex=sp1, sharey=sp1)

energy = np.linspace(0., 2000, num=int(1e3))
sp1.plot(energy, test_obj._delta_e(energy, 1), label=r'$\delta_e$', color='r', linewidth=linewid)
sp1.plot(energy, test_obj._delta_r(energy, 1), label=r'$\delta_r$', color='g', linewidth=linewid)
sp1.plot(energy, test_obj._delta_ts(energy, 1), label=r'$\delta_{ts}$', color='b', linewidth=linewid)
sp1.plot(energy, test_obj._delta_ts(energy, 1) + test_obj._delta_r(energy, 1) + test_obj._delta_e(energy, 1), label=r'$\delta_{tot}$', color='k', linewidth=linewid)
sp1.legend()

for costheta in np.linspace(0, 1, 10):
    sp2.plot(energy, test_obj._delta_e(energy, costheta), label=r'$\delta_e$', color='r', linewidth=linewid)
    sp3.plot(energy, test_obj._delta_r(energy, costheta), label=r'$\delta_r$', color='g', linewidth=linewid)
    sp4.plot(energy, test_obj._delta_ts(energy, costheta), label=r'$\delta_{ts}$', color='b', linewidth=linewid)


for sp in [sp1, sp2, sp3, sp4]:
    sp.grid(alpha=alpha)
plt.suptitle('Furman Pivi test: SEY components', fontsize=30)

# Tests for multiple nn
energy = np.linspace(0.001, E_0_single, num=int(1e5))
E_0 = np.array([E_0_single] * int(1e5))
fig2, axarr = plt.subplots(2, 10, sharex=True, figsize=(1.8 * 12, 12), facecolor='w')
plt.suptitle('Furman-Pivi tests: True secondary energy distributions', fontsize=30)
all_draws = np.array([])
weights_P_n_ts = np.array([])
for kk in np.arange(1, 10.1, 1):
    nn = np.repeat(kk, 1e5)
    prob_density_ts, P_n_ts = test_obj.true_sec_energy_PDF(delta_ts=delta_ts, nn=kk, E_0=E_0_single, energy=energy)
    weights_P_n_ts = np.concatenate([weights_P_n_ts, np.array([P_n_ts] * int(len(energy)))])
    cdf, _ = test_obj.true_sec_energy_CDF(delta_ts=delta_ts, nn=kk, E_0=E_0[0], energy=energy)
    axarr[0, int(kk - 1)].plot(energy, prob_density_ts, label='PDF of true secondary electrons', linewidth=linewid)
    axarr[0, int(kk - 1)].set_title(r'$f_{%i,ts}$' % kk)
    draws = test_obj.get_energy_true_sec(delta_ts=delta_ts, nn=nn, E_0=E_0)
    all_draws = np.concatenate([all_draws, draws])
    axarr[0, int(kk - 1)].hist(draws, density=True, bins=np.arange(0, E_0_single + 1, E_0_single / 100.))
    axarr[0, int(kk - 1)].text(E_0_single / 4., ax4.get_ylim()[1] * 0.3, 'Average: ' + str(np.mean(draws)))
    axarr[0, int(kk - 1)].grid(alpha=alpha)
    axarr[1, int(kk - 1)].plot(energy, cdf, label='CDF', linewidth=linewid)
    axarr[1, int(kk - 1)].grid(alpha=alpha)

plt.figure(10000)
plt.hist(all_draws, density=True, bins=np.arange(0, E_0_single + 1, E_0_single / 100.), weights=weights_P_n_ts)
prob_density_ts = test_obj.average_true_sec_energy_PDF(delta_ts=delta_ts, E_0=E_0_single, energy=energy)
plt.plot(energy, prob_density_ts, 'k', label='PDF of true secondary electrons (average)', linewidth=linewid)
# # Tests for multiple nn
# delta_ts = np.repeat(1.8, 1e3)
# energy = np.linspace(0.001, E_0_single, num=int(1e3))
# energyBig = np.linspace(0.001, E_0_single, num=int(1e5))
# E_0 = np.array([E_0_single] * int(1e3))
# fig, axarr = plt.subplots(2, 10, sharex=True, figsize=(1.8 * 12, 12), facecolor='w')
# nn = np.arange(1, 10.1, 1)
# nn = np.repeat(nn, 10)
# for kk in np.arange(1, 10.1, 1):
#     prob_density_ts, _ = test_obj.true_sec_energy_PDF(delta_ts=delta_ts[0], nn=kk, E_0=E_0[0])
#     axarr[0, int(kk - 1)].plot(energyBig, prob_density_ts, label='PDF of true secondary electrons', linewidth=linewid)
#     axarr[0, int(kk - 1)].set_title(r'$f_{%i,ts}$' % kk)
#     area = scipy.integrate.simps(prob_density_ts, energyBig)
#     area = round(area, round_to_digits)
#     axarr[0, int(kk - 1)].text(50, axarr[0, int(kk - 1)].get_ylim()[1] / 2, 'Area = ' + str(area), fontsize=18)
#     CDF = test_obj.true_sec_energy_CDF(delta_ts=delta_ts[0], nn=kk, E_0=E_0[0])
#     axarr[1, int(kk - 1)].plot(energyBig, CDF, label='CDF', linewidth=linewid)
#     axarr[1, int(kk - 1)].legend()
#
# for kk in np.arange(1, 10.1, 1):
#     nn = np.repeat(kk, 1e3)
#     draw = test_obj.get_energy_true_sec(delta_ts, nn, E_0, choice='poisson', M=10)
#     axarr[0, int(kk - 1)].hist(draw, density=True, bins=np.arange(0, E_0_single + 1, 5))
# fig.subplots_adjust(left=0.05, right=0.95)
plt.show()

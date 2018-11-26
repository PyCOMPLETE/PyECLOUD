import numpy as np
import matplotlib.pyplot as plt
import sec_emission_model_furman_pivi as fp
import mystyle as ms
import scipy

plt.close('all')
ms.mystyle(12)
linewid = 2

test_obj = fp.SEY_model_FP_Cu()  # 276.8, 1.8848)

qq = 0.5  # From FP paper
sigma_e = 2.
E_0 = np.array([300] * int(1e3))
energy = np.linspace(0.001, 300, num=int(1e3))

alpha = 0.3
round_to_digits = 4

plt.close('all')
# Subplots
fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, sharex=True, figsize=(1.8 * 12, 12), facecolor='w')
plt.suptitle('Furman-Pivi tests', fontsize=25)

axlist = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
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
ax1.text(150, ax1.get_ylim()[1] / 2, 'Area = ' + str(area), fontsize=18)
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
ax2.text(150, ax2.get_ylim()[1] / 2, 'Area = ' + str(area), fontsize=18)
ax2.hist(test_obj.get_energy_rediffused(E_0), density=True, bins=20)
ax2.legend()
ax6.legend()

# True secondary
delta_ts = np.repeat(1.8, len(E_0))
prob_density_ts = test_obj.average_true_sec_energy_PDF(delta_ts=1.8, E_0=300, energy=energy)
ax3.plot(energy, prob_density_ts, label='PDF of true secondary electrons', linewidth=linewid)
ax3.set_title('Average true secondary energy distribution')
area = scipy.integrate.simps(prob_density_ts, energy)
area = round(area, round_to_digits)
ax3.text(150, ax3.get_ylim()[1] / 2, 'Area = ' + str(area), fontsize=18)
CDF = test_obj.average_true_sec_energy_CDF(delta_ts=1.8, E_0=300, energy=energy)
ax7.plot(energy, CDF, label='CDF', linewidth=linewid)
# ax7.plot(CDF, energy, label='Inverse CDF')
ax3.hist(test_obj.get_energy_average_true_sec(delta_ts=delta_ts, E_0=E_0, energy=energy), density=True, bins=60)
ax7.legend()

nn = 4
prob_density_ts, _ = test_obj.true_sec_energy_PDF(delta_ts=1.8, nn=nn, E_0=300, energy=energy)
ax4.plot(energy, prob_density_ts, label='PDF of true secondary electrons', linewidth=linewid)
ax4.set_title(r'True secondary energy distribution, $f_{%i,ts}$' % nn)
area = scipy.integrate.simps(prob_density_ts, energy)
area = round(area, round_to_digits)
ax4.text(150, ax4.get_ylim()[1] / 2, 'Area = ' + str(area), fontsize=18)
CDF = test_obj.true_sec_energy_CDF(delta_ts=1.8, nn=nn, E_0=300, energy=energy)
ax8.plot(energy, CDF, label='CDF', linewidth=linewid)
ax4.hist(test_obj.get_energy_true_sec(delta_ts=delta_ts, nn=np.repeat(nn, len(E_0)), E_0=E_0), density=True, bins=30)
ax8.legend()

fig.subplots_adjust(left=0.05, right=0.95)

# SEY components
plt.figure(2, figsize=(9 * 1.2, 9), facecolor='w')
energy = np.linspace(0., 1000, num=int(1e3))
plt.plot(energy, test_obj._delta_e(energy, 1), label=r'$\delta_e$', color='r', linewidth=linewid)
plt.plot(energy, test_obj._delta_r(energy, 1), label=r'$\delta_r$', color='g', linewidth=linewid)
plt.plot(energy, test_obj._delta_ts(energy, 1), label=r'$\delta_{ts}$', color='b', linewidth=linewid)
plt.plot(energy, test_obj._delta_ts(energy, 1) + test_obj._delta_r(energy, 0) + test_obj._delta_e(energy, 0), label=r'$\delta_{tot}$', color='k', linewidth=linewid)
plt.legend()
plt.grid(alpha=alpha)
plt.title('SEY components', fontsize=25)

# # Tests for multiple nn
# delta_ts = np.repeat(1.8, 1e3)
# energy = np.linspace(0.001, 300, num=int(1e3))
# energyBig = np.linspace(0.001, 300, num=int(1e5))
# E_0 = np.array([300] * int(1e3))
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
#     axarr[0, int(kk - 1)].hist(draw, density=True, bins=np.arange(0, 301, 5))
# fig.subplots_adjust(left=0.05, right=0.95)
plt.show()

import numpy as np
import matplotlib.pyplot as plt
import sec_emission_model_furman_pivi as fp
import mystyle as ms
import scipy

plt.close('all')
ms.mystyle(12)

test_obj = fp.SEY_model_FP_Cu()  # 276.8, 1.8848)

qq = 0.5  # From FP paper
sigma_e = 2.
E_0 = np.array([300] * int(1e5))
energy = np.linspace(0.001, 300, num=int(1e5))

alpha = 0.3
round_to_digits = 4

plt.close('all')
# Subplots
fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, sharex=True, figsize=(1.8 * 12, 12))
plt.suptitle('Furman-Pivi tests', fontsize=25)

axlist = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
for ax in axlist:
    ax.grid(alpha=alpha)
    ax.set_xlabel('Energy [eV]')

# Backscattered
prob_density_e = test_obj.backscattered_energy_PDF(energy, E_0)
ax1.plot(energy, prob_density_e, label='PDF')
ax5.plot(energy, test_obj.backscattered_energy_CDF(energy, E_0), label='CDF')
ax1.set_title('Backscattered energy distribution')
area = scipy.integrate.simps(prob_density_e, energy)
area = round(area, round_to_digits)
ax1.text(150, ax1.get_ylim()[1] / 2, 'Area = ' + str(area), fontsize=18)
ax1.legend()
ax5.legend()
ax1.hist(test_obj.get_energy_backscattered(E_0), density=True)

# Rediffused
prob_density_r = test_obj.rediffused_energy_PDF(energy=energy, E_0=E_0)
ax2.plot(energy, prob_density_r, label='PDF')
ax6.plot(energy, test_obj.rediffused_energy_CDF(energy=energy, E_0=E_0), label='CDF')
ax2.set_title('Rediffused energy distribution')
area = scipy.integrate.simps(prob_density_r, energy)
area = round(area, round_to_digits)
ax2.text(150, ax2.get_ylim()[1] / 2, 'Area = ' + str(area), fontsize=18)
ax2.hist(test_obj.get_energy_rediffused(E_0), density=True)
ax2.legend()
ax6.legend()

# True secondary
delta_ts = 1.8
prob_density_ts = test_obj.average_true_sec_energy_PDF(delta_ts=delta_ts, E_0=E_0, choice='binomial')
ax3.plot(energy, prob_density_ts, label='PDF of true secondary electrons')
ax3.set_title('Average true secondary energy distribution')
area = scipy.integrate.simps(prob_density_ts, energy)
area = round(area, round_to_digits)
ax3.text(150, ax3.get_ylim()[1] / 2, 'Area = ' + str(area), fontsize=18)
CDF = test_obj.average_true_sec_energy_CDF(delta_ts=delta_ts, E_0=E_0)
ax7.plot(energy, CDF, label='CDF')
# ax7.plot(CDF, energy, label='Inverse CDF')
ax3.hist(test_obj.get_energy_true_sec(delta_ts=delta_ts, E_0=E_0), density=True)
ax7.legend()

nn = 2
prob_density_ts, _ = test_obj.true_sec_energy_PDF(delta_ts=delta_ts, nn=nn, E_0=E_0, choice='binomial')
ax4.plot(energy, prob_density_ts, label='PDF of true secondary electrons')
ax4.set_title(r'True secondary energy distribution, $f_{%i,ts}$' % nn)
area = scipy.integrate.simps(prob_density_ts, energy)
area = round(area, round_to_digits)
ax4.text(150, ax4.get_ylim()[1] / 2, 'Area = ' + str(area), fontsize=18)

fig.subplots_adjust(left=0.05, right=0.95)

plt.show()

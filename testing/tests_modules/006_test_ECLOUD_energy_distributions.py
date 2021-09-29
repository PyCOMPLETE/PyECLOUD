import sys
if '../..' not in sys.path:
    sys.path.append('../..')
import numpy as np
import matplotlib.pyplot as plt
import sec_emission_model_ECLOUD as ECL
import mystyle as ms
from impact_management_class import impact_management
from geom_impact_ellip import ellip_cham_geom_object
import scipy


plt.close('all')
ms.mystyle(12)
linewid = 2

me = 9.10938356e-31

switch = 2

sey_mod = ECL.SEY_model_ECLOUD(Emax=332., del_max=1.8848, R0=0.7, E_th=35., mufit=1.6636, secondary_angle_distribution='cosine_3D',
                               sigmafit=1.0828, switch_no_increase_energy=switch, thresh_low_energy=1.)

chamb = ellip_cham_geom_object(1., 1., flag_verbose_file=False)
impact_management_object = impact_management(chamb=chamb, sey_mod=sey_mod, Dx_hist=.1, scrub_en_th=25.,
                                             Nbin_En_hist=100, En_hist_max=3000, flag_seg=False,
                                             cos_angle_width=0.05, flag_cos_angle_hist=True)


cos_theta_test = np.linspace(0., 1., 10)
E_0_single = 35.
E_impact_eV_test = E_0_single
n_rep = 100000
alpha = 0.9

extract_ene_hist = impact_management_object.extract_energy_distributions(n_rep, E_impact_eV_test, cos_theta_test, mass=me, Nbin_extract_ene=500, factor_ene_dist_max=1.2)

plt.close('all')
ms.mystyle_arial()

fig1 = plt.figure(1, figsize=(2 * 10, 10))
fig1.set_facecolor('w')
sp1 = fig1.add_subplot(1, 2, 1)
sp2 = fig1.add_subplot(1, 2, 2)


for i_ct, ct in enumerate(cos_theta_test):
    areats = scipy.integrate.simps(extract_ene_hist['true'][:, i_ct], extract_ene_hist['emit_ene_g_hist'])
    areae = scipy.integrate.simps(extract_ene_hist['elast'][:, i_ct], extract_ene_hist['emit_ene_g_hist'])
    thiscol = ms.colorprog(i_ct, len(cos_theta_test))
    label = 'costheta=%.2f' % ct
    sp1.plot(extract_ene_hist['emit_ene_g_hist'], extract_ene_hist['true'][:, i_ct] / areats, color=thiscol, label=label, alpha=alpha, linewidth=linewid, marker='o')
    sp2.plot(extract_ene_hist['emit_ene_g_hist'], extract_ene_hist['elast'][:, i_ct] / areae, color=thiscol, label=label, alpha=alpha, linewidth=linewid, marker='o')

linewid = 3
sp2.plot(0, 0, 'k', label='Model PDF', linewidth=linewid)
sp2.legend(loc='best', prop={'size': 14})
sz = 24
sp1.set_ylabel('True secondaries', fontsize=sz)
sp2.set_ylabel('Elastic', fontsize=sz)

# Compare with model
E_0 = np.array([E_0_single] * int(1e5))
energy = np.linspace(0.001, E_0_single, num=int(1e5))
sigmafit = 1.0828
mufit = 1.6636
hilleret_energy = 1. / (energy * sigmafit * np.sqrt(2 * np.pi)) * np.exp(-(np.log(energy) - mufit)**2 / (2 * sigmafit**2))
area = scipy.integrate.simps(hilleret_energy, energy)
sp1.plot(energy, hilleret_energy / area, 'k', linewidth=linewid)
sp1.set_title('switch_no_increase_energy=%d' % switch)

for sp in [sp1, sp2]:
    sp.grid(alpha=.5)
    sp.set_xlabel('Electron energy [eV]', fontsize=sz)

plt.subplots_adjust(right=0.99, left=.06)

plt.suptitle('Energy distribution extraction tests: ECLOUD model', fontsize=30)

plt.show()

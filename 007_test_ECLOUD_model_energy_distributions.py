import numpy as np
import matplotlib.pyplot as plt
import sec_emission_model_ECLOUD as ECL
import sec_emission_model_furman_pivi as fp
import mystyle as ms
from scipy.constants import e as qe
from impact_management_class import impact_management
from geom_impact_ellip import ellip_cham_geom_object


plt.close('all')
ms.mystyle(12)
linewid = 2

me = 9.10938356e-31

furman_pivi_surface_LHC = {'M': 2,
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

sey_mod = fp.SEY_model_furman_pivi(E_th=35., sigmafit=1.0828, mufit=1.6636, secondary_angle_distribution='cosine_3D',
                                   switch_no_increase_energy=0, thresh_low_energy=-1,
                                   furman_pivi_surface=furman_pivi_surface_LHC)
sey_mod = ECL.SEY_model_ECLOUD(Emax=332., del_max=1.8848, R0=0.7, E_th=35., mufit=1.6636, secondary_angle_distribution='cosine_3D',
                               sigmafit=1.0828, switch_no_increase_energy=0, thresh_low_energy=-1)


# def extract_energy_distributions(n_rep, E_impact_eV_test, cos_theta_test, charge, mass):
#     dists = {}
#     for etype in sey_mod.event_types.keys():
#         etype_name = sey_mod.event_types[etype]
#         dists[etype_name] = []
#     print('Extracting energy distributions...')
#     for i_ct, ct in enumerate(cos_theta_test):
#         print('%d/%d' % (i_ct + 1, len(cos_theta_test)))
#         Ene = E_impact_eV_test
#         nel_impact = np.ones(n_rep)
#         # Assuming normal is along x
#         v_mod = np.sqrt(2 * Ene * qe / mass) * np.ones_like(nel_impact)
#         vx = v_mod * ct
#         vy = v_mod * np.sqrt(1 - ct * ct)
#
#         nel_emit_tot_events, event_type, event_info,\
#             nel_replace, x_replace, y_replace, z_replace, vx_replace, vy_replace, vz_replace, i_seg_replace,\
#             nel_new_MPs, x_new_MPs, y_new_MPs, z_new_MPs, vx_new_MPs, vy_new_MPs, vz_new_MPs, i_seg_new_MPs =\
#             sey_mod.impacts_on_surface(
#                 mass=mass, nel_impact=nel_impact, x_impact=nel_impact * 0, y_impact=nel_impact * 0, z_impact=nel_impact * 0,
#                 vx_impact=vx * np.ones_like(nel_impact),
#                 vy_impact=vy * np.ones_like(nel_impact),
#                 vz_impact=nel_impact * 0,
#                 Norm_x=np.ones_like(nel_impact), Norm_y=np.zeros_like(nel_impact),
#                 i_found=np.int_(np.ones_like(nel_impact)),
#                 v_impact_n=vx * np.ones_like(nel_impact),
#                 E_impact_eV=Ene * np.ones_like(nel_impact),
#                 costheta_impact=ct * np.ones_like(nel_impact),
#                 nel_mp_th=1,
#                 flag_seg=True)
#
#         v_replace_mod = np.sqrt(vx_replace**2 + vy_replace**2 + vz_replace**2)
#         E_replace_eV = 0.5 * mass / qe * v_replace_mod * v_replace_mod
#
#         v_new_MPs_mod = np.sqrt(vx_new_MPs**2 + vy_new_MPs**2 + vz_new_MPs**2)
#         E_new_MPs_eV = 0.5 * mass / qe * v_new_MPs_mod * v_new_MPs_mod
#
#         E_all_MPs_eV = np.concatenate([E_replace_eV, E_new_MPs_eV])
#         extended_nel_emit_tot_events = event_info['extended_nel_emit_tot_events']
#         E_all_MPs_eV = E_all_MPs_eV * extended_nel_emit_tot_events
#
#         extended_event_type = event_info['extended_event_type']
#         for etype in sey_mod.event_types.keys():
#             etype_name = sey_mod.event_types[etype]
#             dists[etype_name].append(E_all_MPs_eV[extended_event_type == etype])
#
#     print('Done extracting energy distributions.')
#
#     return dists


chamb = ellip_cham_geom_object(1., 1., flag_verbose_file=False)
impact_management_object = impact_management(chamb=chamb, sey_mod=sey_mod, Dx_hist=.1, scrub_en_th=25.,
                                             Nbin_En_hist=100, En_hist_max=3000, flag_seg=False,
                                             cos_angle_width=0.05, flag_cos_angle_hist=True)


cos_theta_test = np.linspace(1., 1., 1)
E_0_single = 100
E_impact_eV_test = np.array([E_0_single] * int(1e5))
n_rep = 100000
alpha = 0.9

dists = impact_management_object.extract_energy_distributions(n_rep, E_impact_eV_test, cos_theta_test, charge=qe, mass=me)

plt.close('all')
ms.mystyle_arial()

fig1 = plt.figure(1, figsize=(2 * 8, 8))
fig1.set_facecolor('w')
sp1 = fig1.add_subplot(1, 2, 1)
sp2 = fig1.add_subplot(1, 2, 2)


for i_ct, ct in enumerate(cos_theta_test):
    thiscol = ms.colorprog(i_ct, len(cos_theta_test))
    label = 'costheta=%.2f' % ct
    sp1.hist(dists['true'][i_ct], bins=60, color=thiscol, label=label, alpha=alpha, density=True)
    sp2.hist(dists['elast'][i_ct], bins=30, color=thiscol, label=label, alpha=alpha, density=True)

linewid = 3
sp2.plot(0, 0, 'k', label='Model PDF', linewidth=linewid)
sp2.legend(loc='best', prop={'size': 14})
sz = 24
sp1.set_ylabel('True secondaries', fontsize=sz)
sp2.set_ylabel('Elastic', fontsize=sz)

# Compare with model
E_0 = np.array([E_0_single] * int(1e5))
energy = np.linspace(0.001, E_0_single, num=int(1e5))


for sp in [sp1, sp2]:
    sp.grid('on')
    sp.set_xlabel('Electron energy [eV]')

plt.subplots_adjust(right=0.99, left=.06)

plt.suptitle('Energy distribution extraction tests: ECLOUD model', fontsize=30)

plt.show()

import numpy as np
import matplotlib.pyplot as plt
import sec_emission_model_ECLOUD as ECL
import sec_emission_model_furman_pivi as fp
import mystyle as ms
from scipy.constants import e as qe

plt.close('all')
ms.mystyle(12)
linewid = 2

me = 9.10938356e-31
sey_mod = fp.SEY_model_FP_Cu(E_th=35., secondary_angle_distribution='cosine_3D',
                             switch_no_increase_energy=0, thresh_low_energy=-1)  # 276.8, 1.8848)
# sey_mod = ECL.SEY_model_ECLOUD(Emax=332., del_max=1.8848, R0=0.7, E_th=35., mufit=1.6636, secondary_angle_distribution='cosine_3D',
#                                sigmafit=1.0828, switch_no_increase_energy=0, thresh_low_energy=-1)


def extract_sey_curves(n_rep, E_impact_eV_test, cos_theta_test, charge, mass):

    # del_true_mat = np.zeros((len(cos_theta_test), len(E_impact_eV_test)))
    # del_elast_mat = np.zeros((len(cos_theta_test), len(E_impact_eV_test)))
    # del_rediff_mat = np.zeros((len(cos_theta_test), len(E_impact_eV_test)))
    # del_absorb_mat = np.zeros((len(cos_theta_test), len(E_impact_eV_test)))
    # del_tot_mat = np.zeros((len(cos_theta_test), len(E_impact_eV_test)))

    deltas = {}
    for etype in sey_mod.event_types.keys():
        etype_name = sey_mod.event_types[etype]
        deltas[etype_name] = np.zeros((len(cos_theta_test), len(E_impact_eV_test)))
    print('Extracting SEY curves...')
    for i_ct, ct in enumerate(cos_theta_test):
        print('%d/%d' % (i_ct + 1, len(cos_theta_test)))
        for i_ene, Ene in enumerate(E_impact_eV_test):

            # nel_emit, flag_elast, flag_truesec = sey_mod.SEY_process(nel_impact=np.ones(n_rep),
            #                 E_impact_eV=Ene*np.ones(n_rep), costheta_impact=np.ones(n_rep)*ct, i_impact=np.array(n_rep*[0]))
            nel_impact = np.ones(n_rep)
            # Assuming normal is along x
            v_mod = np.sqrt(2 * Ene * qe / mass) * np.ones_like(nel_impact)
            vx = v_mod * ct
            vy = v_mod * np.sqrt(1 - ct * ct)

            nel_emit_tot_events, event_type, event_info,\
                nel_replace, x_replace, y_replace, z_replace, vx_replace, vy_replace, vz_replace, i_seg_replace,\
                nel_new_MPs, x_new_MPs, y_new_MPs, z_new_MPs, vx_new_MPs, vy_new_MPs, vz_new_MPs, i_seg_new_MPs =\
                sey_mod.impacts_on_surface(
                    mass=mass, nel_impact=nel_impact, x_impact=nel_impact * 0, y_impact=nel_impact * 0, z_impact=nel_impact * 0,
                    vx_impact=vx * np.ones_like(nel_impact),
                    vy_impact=vy * np.ones_like(nel_impact),
                    vz_impact=nel_impact * 0,
                    Norm_x=np.ones_like(nel_impact), Norm_y=np.zeros_like(nel_impact),
                    i_found=np.int_(np.ones_like(nel_impact)),
                    v_impact_n=vx * np.ones_like(nel_impact),
                    E_impact_eV=Ene * np.ones_like(nel_impact),
                    costheta_impact=ct * np.ones_like(nel_impact),
                    nel_mp_th=1,
                    flag_seg=True)
            # del_tot = (np.sum(nel_replace) + np.sum(nel_new_MPs)) / float(n_rep)
            # del_tot_mat[i_ct, i_ene] = del_tot
            # # del_true_mat[i_ct, i_ene] = del_tot * (float(np.sum(event_type == 1)) + float(np.sum(event_type == 3))) / float(n_rep)
            # del_true_mat[i_ct, i_ene] = (np.sum(nel_replace) + np.sum(nel_new_MPs) - np.sum(event_type == 0) - np.sum(event_type == 2)) / float(n_rep)
            # del_tot = 1.
            # del_elast_mat[i_ct, i_ene] = del_tot * float(np.sum(event_type == 0)) / float(n_rep)
            # del_rediff_mat[i_ct, i_ene] = del_tot * float(np.sum(event_type == 2)) / float(n_rep)
            # del_absorb_mat[i_ct, i_ene] = del_tot * float(np.sum(event_type == 3)) / float(n_rep)

            for etype in sey_mod.event_types.keys():
                etype_name = sey_mod.event_types[etype]
                thisdelta = deltas[etype_name]
                thisdelta[i_ct, i_ene] = np.sum(nel_emit_tot_events[event_type == etype]) / np.sum(nel_impact)
                deltas[etype_name] = thisdelta

    print('Done extracting SEY curves.')

    return deltas


cos_theta_test = np.linspace(0, 1., 10)
E_impact_eV_test = np.array(list(np.arange(0, 499., 1.)) + list(np.arange(500., 2000, 5.)))
n_rep = 10000

deltas = extract_sey_curves(n_rep, E_impact_eV_test, cos_theta_test, charge=qe, mass=me)
del_true_mat = deltas['true']
del_elast_mat = deltas['elast']
del_rediff_mat = deltas['rediff']
del_absorb_mat = deltas['absorb']

plt.close('all')
ms.mystyle_arial()

fig1 = plt.figure(1, figsize=(3 * 8, 2 * 8))
fig1.set_facecolor('w')
sp1 = fig1.add_subplot(2, 3, 1)
sp2 = fig1.add_subplot(2, 3, 2, sharex=sp1)
sp3 = fig1.add_subplot(2, 3, 3, sharex=sp1)
sp4 = fig1.add_subplot(2, 3, 4, sharex=sp1)
sp5 = fig1.add_subplot(2, 3, 5, sharex=sp1)
sp6 = fig1.add_subplot(2, 3, 6, sharex=sp1)

for i_ct, ct in enumerate(cos_theta_test):
    thiscol = ms.colorprog(i_ct, len(cos_theta_test))
    label = 'costheta=%.2f' % ct
    sp1.plot(E_impact_eV_test, del_true_mat[i_ct, :], color=thiscol, label=label)
    sp2.plot(E_impact_eV_test, del_elast_mat[i_ct, :], color=thiscol, label=label)
    sp3.plot(E_impact_eV_test, del_rediff_mat[i_ct, :], color=thiscol, label=label)
    sp4.plot(E_impact_eV_test, del_absorb_mat[i_ct, :], color=thiscol, label=label)
    sp5.plot(E_impact_eV_test, del_true_mat[i_ct, :] + del_rediff_mat[i_ct, :] + del_elast_mat[i_ct, :], color=thiscol, label=label)
    sp6.plot(E_impact_eV_test, del_true_mat[i_ct, :] + del_elast_mat[i_ct, :], color=thiscol, label=label)

sp3.plot(0, 0, 'white', label='Model')
sp3.legend(loc='best', prop={'size': 14})
sp1.set_ylabel('Delta true')
sp2.set_ylabel('Delta elast')
sp3.set_ylabel('Delta rediff')
sp4.set_ylabel('Delta absorb')
sp5.set_ylabel('Delta total')
sp6.set_ylabel(r'$\delta_{ts} + \delta_{e}$')

for sp in [sp1, sp2, sp3, sp4, sp5, sp6]:
    sp.grid('on')
    sp.set_xlabel('Electron energy [eV]')

plt.subplots_adjust(right=0.99, left=.05)


test_obj = fp.SEY_model_FP_Cu()  # 276.8, 1.8848)

energy = np.linspace(0., 2000, num=int(1e3))

for costheta in np.linspace(0, 1, 10):
    delta_ts_vec = test_obj._delta_ts(energy, costheta)
    delta_e_vec = test_obj._delta_e(energy, costheta)
    delta_r_vec = test_obj._delta_r(energy, costheta)

    sp2.plot(energy, delta_e_vec, color='k', linewidth=linewid)
    sp3.plot(energy, delta_r_vec, color='k', linewidth=linewid)
    sp1.plot(energy, delta_ts_vec, color='k', linewidth=linewid)
    sp5.plot(energy, delta_r_vec + delta_ts_vec + delta_e_vec, color='k', linewidth=linewid)
    sp6.plot(energy, delta_ts_vec + delta_e_vec, color='k', linewidth=linewid)

plt.suptitle('SEY extraction tests: Furman-Pivi model', fontsize=30)

plt.show()

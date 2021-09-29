import sys
if '../..' not in sys.path:
    sys.path.append('../..')
import numpy as np
import matplotlib.pyplot as plt
import sec_emission_model_ECLOUD as ECL
import mystyle as ms
from scipy.constants import e as qe

plt.close('all')
ms.mystyle(12)
linewid = 2
del_max = 1.8
me = 9.10938356e-31


def del_elas_ECLOUD(energy, cos_theta, R_0=0.7, E_max=332., E_0=150.):
    del_elas = R_0 * ((np.sqrt(energy) - np.sqrt(energy + E_0)) / (np.sqrt(energy) + np.sqrt(energy + E_0)))**2
    return del_elas


def del_true_ECLOUD(energy, cos_theta, del_max, s=1.35, E_max=332.):
    E_max_tilde = E_max * (1. + 0.7 * (1. - costheta))
    x = energy / E_max_tilde
    del_true = del_max * s * x / (s - 1 + x**s)
    return del_true * np.exp(0.5 * (1. - costheta))


switch = 0
sey_mod = ECL.SEY_model_ECLOUD(Emax=332., del_max=del_max, R0=0.7, E_th=35., mufit=1.6636, secondary_angle_distribution='cosine_3D',
                               sigmafit=1.0828, switch_no_increase_energy=switch, thresh_low_energy=1.,
                               flag_costheta_delta_scale=True, flag_costheta_Emax_shift=True)


def extract_sey_curves(n_rep, E_impact_eV_test, cos_theta_test, charge, mass):

    deltas = {}
    for etype in list(sey_mod.event_types.keys()):
        etype_name = sey_mod.event_types[etype]
        deltas[etype_name] = np.zeros((len(cos_theta_test), len(E_impact_eV_test)))
    print('Extracting SEY curves...')
    for i_ct, ct in enumerate(cos_theta_test):
        print(('%d/%d' % (i_ct + 1, len(cos_theta_test))))
        for i_ene, Ene in enumerate(E_impact_eV_test):

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

            for etype in list(sey_mod.event_types.keys()):
                etype_name = sey_mod.event_types[etype]
                thisdelta = deltas[etype_name]
                thisdelta[i_ct, i_ene] = np.sum(
                    nel_emit_tot_events[event_type == etype]) / np.sum(nel_impact)
                deltas[etype_name] = thisdelta

    print('Done extracting SEY curves.')

    return deltas


cos_theta_test = np.linspace(0, 1., 10)
E_impact_eV_test = np.array(list(np.arange(0, 499., 5.)) + list(np.arange(500., 2000, 25.)))
n_rep = int(1e3)

deltas = extract_sey_curves(n_rep, E_impact_eV_test, cos_theta_test, charge=qe, mass=me)
del_true_mat = deltas['true']
del_elast_mat = deltas['elast']

plt.close('all')
ms.mystyle_arial()

fig1 = plt.figure(1, figsize=(3 * 8, 2 * 8))
fig1.set_facecolor('w')
sp1 = fig1.add_subplot(1, 3, 1)
sp2 = fig1.add_subplot(1, 3, 2, sharex=sp1)
sp5 = fig1.add_subplot(1, 3, 3, sharex=sp1)

for i_ct, ct in enumerate(cos_theta_test):
    thiscol = ms.colorprog(i_ct, len(cos_theta_test))
    label = 'costheta=%.2f' % ct
    sp1.plot(E_impact_eV_test, del_true_mat[i_ct, :], color=thiscol, label=label, linewidth=linewid)
    sp2.plot(E_impact_eV_test, del_elast_mat[i_ct, :], color=thiscol, label=label, linewidth=linewid)
    sp5.plot(E_impact_eV_test, del_true_mat[i_ct, :] + del_elast_mat[i_ct, :], color=thiscol, label=label, linewidth=linewid)

sp1.set_ylabel('Delta true')
sp2.set_ylabel('Delta elast')
sp5.set_ylabel('Delta total')

for sp in [sp1, sp2, sp5]:
    sp.grid('on')
    sp.set_xlabel('Electron energy [eV]')

plt.subplots_adjust(right=0.99, left=.05)

test_obj = sey_mod

energy = np.linspace(0., 2000, num=int(1e3))

for costheta in np.linspace(0, 1., 10):
    delta_ts_vec = del_true_ECLOUD(energy, costheta, del_max=del_max)
    delta_e_vec = del_elas_ECLOUD(energy, costheta)

    sp2.plot(energy, delta_e_vec, color='k', linewidth=linewid, linestyle='--')
    sp1.plot(energy, delta_ts_vec, color='k', linewidth=linewid, linestyle='--')
    sp5.plot(energy, delta_ts_vec + delta_e_vec, color='k', linewidth=linewid, linestyle='--')

sp2.legend(loc='best', prop={'size': 14})

plt.suptitle('SEY extraction tests: ECLOUD model', fontsize=30)

plt.show()

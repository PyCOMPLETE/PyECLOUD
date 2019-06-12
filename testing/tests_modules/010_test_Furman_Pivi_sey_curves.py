import sys
if '../..' not in sys.path:
    sys.path.append('../..')
import numpy as np
import matplotlib.pyplot as plt
import sec_emission_model_furman_pivi as fp
import mystyle as ms
from scipy.constants import e as qe

plt.close('all')
ms.mystyle(12)
linewid = 2

me = 9.10938356e-31


def del_elas_ECLOUD(energy, R_0=0.7, E_max=332., E_0=150.):
    del_elas = R_0 * ((np.sqrt(energy) - np.sqrt(energy + E_0)) / (np.sqrt(energy) + np.sqrt(energy + E_0)))**2
    return del_elas


def del_true_ECLOUD(energy, del_max, s=1.35, E_max=332., costheta=1.):
    angular_factor = np.exp(0.5 * (1. - costheta))
    E_max_tilde = E_max * (1. + 0.7 * (1. - costheta))
    x = energy / E_max_tilde
    del_true = del_max * s * x / (s - 1 + x**s)
    return del_true * angular_factor


furman_pivi_surface_tweak = {
    'use_ECLOUD_theta0_dependence': True,
    'use_ECLOUD_energy': False,
    'conserve_energy': False,
    'exclude_rediffused': True,
    'choice': 'poisson',
    'M_cut': 10,
    'p_n': np.array([1.21963859, 1.66070543, 1.21935223, 1.09987752, 4.28158656, 1.02052557, 1.0247471, 1.02307995, 29.93491271, 1.02045612]),
    'eps_n': np.array([7.44033631e+00, 2.47339424e+00, 7.45004962e+00, 1.63618903e+01, 4.97986255e-01, 7.96170380e+01, 6.60354258e+01, 7.08053955e+01, 5.64779654e-02, 7.98873331e+01]),
    # Parameters for backscattered electrons
    'p1EInf': 0.002158,  # Changed this
    'p1Ehat': 0.709633,  # Changed this
    'eEHat': 0.,
    'w': 46.028959,  # Changed this
    'p': 0.468907,  # Changed this
    'e1': 0.,  # Changed this
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
    'deltaTSHat': 1.8,
    'eHat0': 332.,
    's': 1.35,
    't1': 0.5, # 0.706340,  # Changed this
    't2': 1., #0.715223,  # Changed this
    't3': 0.7,
    't4': 1.}

furman_pivi_surface_LHC = {
    'use_ECLOUD_theta0_dependence': False,
    'use_ECLOUD_energy': False,
    'conserve_energy': False,
    'exclude_rediffused': False,
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
    't4': 1.}

furman_pivi_surface = {
    'use_ECLOUD_theta0_dependence': False,
    'use_ECLOUD_energy': False,
    'conserve_energy': False,
    'exclude_rediffused': False,
    'choice': 'poisson',
    'M_cut': 10,
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
    'deltaTSHat': 1.6 - 0.22, #1.8848,
    'eHat0': 276.8,
    's': 1.54,
    't1': 0.66,
    't2': 0.8,
    't3': 0.7,
    't4': 1.}

# Scaled py POSINST to del_tot_max = 1.6
furman_pivi_surface_scaled = {
    'use_ECLOUD_theta0_dependence': False,
    'use_ECLOUD_energy': False,
    'conserve_energy': False,
    'exclude_rediffused': False,
    'choice': 'poisson',
    'M_cut': 10,
    'p_n': np.array([2.5, 3.3, 2.5, 2.5, 2.8, 1.3, 1.5, 1.5, 1.5, 1.5]),
    'eps_n': np.array([1.5, 1.75, 1., 3.75, 8.5, 11.5, 2.5, 3., 2.5, 3.]),
    'p1EInf': 0.015294,  # Changed this
    'p1Ehat': 0.382362,  # Changed this
    'eEHat': 0.,
    'w': 60.86,
    'p': 1.,
    'e1': 0.26,
    'e2': 2.,
    'sigmaE': 2.,
    'p1RInf': 0.152945,  # Changed this
    'eR': 0.041,
    'r': 0.104,
    'q': 0.5,
    'r1': 0.26,
    'r2': 2.,
    'deltaTSHat': 1.441353,  # Changed this
    'eHat0': 276.8,
    's': 1.54,
    't1': 0.66,
    't2': 0.8,
    't3': 0.7,
    't4': 1.}

flag_costheta_delta_scale = True
flag_costheta_Emax_shift = True

sey_mod = fp.SEY_model_furman_pivi(E_th=35., sigmafit=1.0828, mufit=1.6636, secondary_angle_distribution='cosine_3D',
                                   switch_no_increase_energy=0, thresh_low_energy=-1,
                                   furman_pivi_surface=furman_pivi_surface, flag_costheta_delta_scale=flag_costheta_delta_scale,
                                   flag_costheta_Emax_shift=flag_costheta_Emax_shift)


def extract_sey_curves(n_rep, E_impact_eV_test, cos_theta_test, charge, mass):

    deltas = {}
    for etype in sey_mod.event_types.keys():
        etype_name = sey_mod.event_types[etype]
        deltas[etype_name] = np.zeros((len(cos_theta_test), len(E_impact_eV_test)))
    print('Extracting SEY curves...')
    for i_ct, ct in enumerate(cos_theta_test):
        print('%d/%d' % (i_ct + 1, len(cos_theta_test)))
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

            for etype in sey_mod.event_types.keys():
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
    sp1.plot(E_impact_eV_test, del_true_mat[i_ct, :], color=thiscol, label=label, linewidth=linewid)
    sp2.plot(E_impact_eV_test, del_elast_mat[i_ct, :], color=thiscol, label=label, linewidth=linewid)
    sp3.plot(E_impact_eV_test, del_rediff_mat[i_ct, :], color=thiscol, label=label, linewidth=linewid)
    sp4.plot(E_impact_eV_test, del_absorb_mat[i_ct, :], color=thiscol, label=label, linewidth=linewid)
    sp5.plot(E_impact_eV_test, del_true_mat[i_ct, :] + del_rediff_mat[i_ct, :] + del_elast_mat[i_ct, :], color=thiscol, label=label, linewidth=linewid)
    sp6.plot(E_impact_eV_test, del_true_mat[i_ct, :] + del_elast_mat[i_ct, :], color=thiscol, label=label, linewidth=linewid)

sp3.plot(0, 0, 'white', label='Model')
sp1.set_ylabel('Delta true')
sp2.set_ylabel('Delta elast')
sp3.set_ylabel('Delta rediff')
sp4.set_ylabel('Delta absorb')
sp5.set_ylabel('Delta total')
sp6.set_ylabel(r'$\delta_{ts} + \delta_{e}$')

for sp in [sp1, sp2, sp3, sp4, sp5, sp6]:
    sp.grid(True)
    sp.set_xlabel('Electron energy [eV]')

plt.subplots_adjust(right=0.99, left=.05)

test_obj = sey_mod

energy = np.linspace(0., 2000, num=int(1e5))

for costheta in np.linspace(0, 1., 10):
    delta_ts_vec = test_obj.delta_ts(energy, costheta)
    delta_e_vec = test_obj.delta_e(energy, costheta)
    delta_r_vec = test_obj.delta_r(energy, costheta)

    sp2.plot(energy, delta_e_vec, color='k', linewidth=linewid)
    sp3.plot(energy, delta_r_vec, color='k', linewidth=linewid)
    sp1.plot(energy, delta_ts_vec, color='k', linewidth=linewid)
    sp5.plot(energy, delta_r_vec + delta_ts_vec + delta_e_vec, color='k', linewidth=linewid)
    sp6.plot(energy, delta_ts_vec + delta_e_vec, color='k', linewidth=linewid)

sp2.plot(energy, del_elas_ECLOUD(energy), '--', color='r', linewidth=linewid, label='ECLOUD model')
for ct in cos_theta_test:
    sp1.plot(energy, del_true_ECLOUD(energy, del_max=test_obj.deltaTSHat, costheta=ct), '--', color='r', linewidth=linewid, label='ECLOUD model')
sp2.legend(loc='best', prop={'size': 14})

plt.suptitle('SEY extraction tests: Furman-Pivi model \nexclude_rediffused=%s' % str(sey_mod.exclude_rediffused), fontsize=30)

plt.show()

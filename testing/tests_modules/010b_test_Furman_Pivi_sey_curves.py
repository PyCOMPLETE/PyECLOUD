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
ms.mystyle(25)
sz = 30


def del_elas_ECLOUD(energy, R_0=0.7, E_max=332., E_0=150.):
    del_elas = R_0 * ((np.sqrt(energy) - np.sqrt(energy + E_0)) / (np.sqrt(energy) + np.sqrt(energy + E_0)))**2
    return del_elas


def del_true_ECLOUD(energy, del_max, s=1.35, E_max=332.):
    x = energy / E_max
    del_true = del_max * s * x / (s - 1 + x**s)
    return del_true


furman_pivi_surface_tweak = {
    'use_modified_sigmaE': False,
    'use_ECLOUD_theta0_dependence': False,
    'use_ECLOUD_energy': False,
    'conserve_energy': True,
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
    'deltaTSHat': 1.8848,
    'eHat0': 332.,
    's': 1.35,
    't1': 0.706340,  # Changed this
    't2': 0.715223,  # Changed this
    't3': 0.7,
    't4': 1.}

furman_pivi_surface_LHC = {
    'use_modified_sigmaE': False,
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
    'use_modified_sigmaE': False,
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
    'deltaTSHat': 1.8848,
    'eHat0': 276.8,
    's': 1.54,
    't1': 0.66,
    't2': 0.8,
    't3': 0.7,
    't4': 1.}

sey_mod = fp.SEY_model_furman_pivi(E_th=35., sigmafit=1.0828, mufit=1.6636, secondary_angle_distribution='cosine_3D',
                                   switch_no_increase_energy=0, thresh_low_energy=-1,
                                   furman_pivi_surface=furman_pivi_surface_LHC)


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


E_ts = 3000.
E_e = 500.
E_r = 10.
cos_theta_test = np.linspace(0, 1., 10)
E_impact_eV_test = np.array(list(np.arange(0, 499., 5.)) + list(np.arange(500., E_ts, 25.)))
E_impact_eV_test_e = np.array(list(np.arange(0, E_e, .5)))
E_impact_eV_test_r = np.array(list(np.arange(0, E_r, .05)))
n_rep = int(1e3)

deltas = extract_sey_curves(n_rep, E_impact_eV_test, cos_theta_test, charge=qe, mass=me)
deltas_e = extract_sey_curves(n_rep, E_impact_eV_test_e, cos_theta_test, charge=qe, mass=me)
deltas_r = extract_sey_curves(n_rep, E_impact_eV_test_r, cos_theta_test, charge=qe, mass=me)

del_true_mat = deltas['true']
del_elast_mat = deltas['elast']
del_rediff_mat = deltas['rediff']
del_absorb_mat = deltas['absorb']

del_elast_mat_2 = deltas_e['elast']
del_rediff_mat_2 = deltas_r['rediff']


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

fig_ts = plt.figure(2, figsize=(12, 9), facecolor='white')
plt.xticks(fontsize=sz - 5)
plt.yticks(fontsize=sz - 5)
plt.subplots_adjust(bottom=.11)
sp_ts = fig_ts.add_subplot(1, 1, 1)
fig_e = plt.figure(3, figsize=(12, 9), facecolor='white')
plt.xticks(fontsize=sz - 5)
plt.yticks(fontsize=sz - 5)
plt.subplots_adjust(bottom=.11)
sp_e = fig_e.add_subplot(1, 1, 1)
fig_r = plt.figure(4, figsize=(12, 9), facecolor='white')
plt.xticks(fontsize=sz - 5)
plt.yticks(fontsize=sz - 5)
plt.subplots_adjust(bottom=.11)
sp_r = fig_r.add_subplot(1, 1, 1)

plt.figure(fig1.number)
plt.subplots_adjust(right=0.99, left=.07)

test_obj = sey_mod

energy = np.linspace(0., E_ts, num=int(1e3))
energy_e = np.linspace(0., E_e, num=int(2e2))
energy_r = np.linspace(0., E_r, num=int(1e3))

for i_ct, ct in enumerate(cos_theta_test):
    thiscol = ms.colorprog(i_ct, len(cos_theta_test))
    label = 'costheta=%.2f' % ct
    sp1.plot(E_impact_eV_test, del_true_mat[i_ct, :], color=thiscol, label=label, linewidth=linewid)
    sp_ts.plot(E_impact_eV_test, del_true_mat[i_ct, :], color=thiscol, label=label, linewidth=linewid)
    sp2.plot(E_impact_eV_test, del_elast_mat[i_ct, :], color=thiscol, label=label, linewidth=linewid)
    sp_e.plot(E_impact_eV_test_e, del_elast_mat_2[i_ct, :], color=thiscol, label=label, linewidth=linewid)
    sp3.plot(E_impact_eV_test, del_rediff_mat[i_ct, :], color=thiscol, label=label, linewidth=linewid)
    sp_r.plot(E_impact_eV_test_r, del_rediff_mat_2[i_ct, :], color=thiscol, label=label, linewidth=linewid)
    sp4.plot(E_impact_eV_test, del_absorb_mat[i_ct, :], color=thiscol, label=label, linewidth=linewid)
    sp5.plot(E_impact_eV_test, del_true_mat[i_ct, :] + del_rediff_mat[i_ct, :] + del_elast_mat[i_ct, :], color=thiscol, label=label, linewidth=linewid)
    sp6.plot(E_impact_eV_test, del_true_mat[i_ct, :] + del_elast_mat[i_ct, :], color=thiscol, label=label, linewidth=linewid)

stl = '--'
for costheta in np.linspace(0, 1., 10):
    delta_ts_vec = test_obj.delta_ts(energy, costheta)
    delta_e_vec = test_obj.delta_e(energy, costheta)
    delta_r_vec = test_obj.delta_r(energy, costheta)

    delta_e_vec_2 = test_obj.delta_e(energy_e, costheta)
    delta_r_vec_2 = test_obj.delta_r(energy_r, costheta)

    sp2.plot(energy, delta_e_vec, color='k', linewidth=linewid, linestyle=stl)
    sp_e.plot(energy_e, delta_e_vec_2, color='k', linewidth=linewid, linestyle=stl)
    sp3.plot(energy, delta_r_vec, color='k', linewidth=linewid, linestyle=stl)
    sp_r.plot(energy_r, delta_r_vec_2, color='k', linewidth=linewid, linestyle=stl)
    sp1.plot(energy, delta_ts_vec, color='k', linewidth=linewid, linestyle=stl)
    sp_ts.plot(energy, delta_ts_vec, color='k', linewidth=linewid, linestyle=stl)
    sp5.plot(energy, delta_r_vec + delta_ts_vec + delta_e_vec, color='k', linewidth=linewid, linestyle=stl)
    sp6.plot(energy, delta_ts_vec + delta_e_vec, color='k', linewidth=linewid, linestyle=stl)

sp3.plot(0, 0, 'white', label='Model')
sp1.set_ylabel(r'$\delta_{ts}$', fontsize=sz)
sp_ts.set_ylabel(r'$\delta_{ts}$', fontsize=sz)
sp2.set_ylabel(r'$\delta_{e}$', fontsize=sz)
sp_e.set_ylabel(r'$\delta_{e}$', fontsize=sz)
sp3.set_ylabel(r'$\delta_{r}$', fontsize=sz)
sp_r.set_ylabel(r'$\delta_{r}$', fontsize=sz)
sp4.set_ylabel(r'$\delta_{absorbed}$', fontsize=sz)
sp5.set_ylabel(r'$\delta_{tot}$', fontsize=sz)
sp6.set_ylabel(r'$\delta_{ts} + \delta_{e}$', fontsize=sz)

for sp in [sp1, sp2, sp3, sp4, sp5, sp6, sp_ts, sp_e, sp_r]:
    sp.grid('on')
    sp.set_xlabel('Electron energy [eV]', fontsize=sz)


# sp2.plot(energy, del_elas_ECLOUD(energy), '--', color='r', linewidth=linewid, label='ECLOUD model \nnormal incidence')
# sp1.plot(energy, del_true_ECLOUD(energy, del_max=furman_pivi_surface_LHC['deltaTSHat']), '--', color='r', linewidth=linewid, label='ECLOUD model')
# sp2.legend(loc='best', prop={'size': 14})

plt.suptitle('SEY extraction tests: Furman-Pivi model \nexclude_rediffused=%s'%str(furman_pivi_surface_LHC['exclude_rediffused']), fontsize=30)

plt.show()

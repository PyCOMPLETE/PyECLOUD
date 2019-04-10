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

furman_pivi_surface_tweak = {'use_ECLOUD_energy': False,
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
                             't4': 1.,
                             }

furman_pivi_surface_LHC = {'use_ECLOUD_energy': False,
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
                           't4': 1.,
                           }

furman_pivi_surface = {'use_ECLOUD_energy': False,
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
                       't4': 1.,
                       }


secondary_angle_distribution = 'cosine_3D'

sey_mod = fp.SEY_model_furman_pivi(E_th=35., sigmafit=1.0828, mufit=1.6636, secondary_angle_distribution=secondary_angle_distribution,
                                   switch_no_increase_energy=0, thresh_low_energy=-1,
                                   furman_pivi_surface=furman_pivi_surface_tweak)


def extract_emission_angles(n_rep, E_impact_eV_test, cos_theta_test, charge, mass):
    """
    Extracts emission angles as generated from sey_mod.impacts_on_surface.
    Returns an array containing n_rep emission angles corresponding to n_rep MPs
    that have been allowed to impact on the surface with an incident angle
    cos_theta_test and impacting energy E_impact_eV_test.
    """
    print('Extracting Emission angles...')
    Ene = E_impact_eV_test
    ct = cos_theta_test

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

    # Computing emission angles
    vx_all = np.concatenate((vx_replace, vx_new_MPs))
    vy_all = np.concatenate((vy_replace, vy_new_MPs))
    vz_all = np.concatenate((vz_replace, vz_new_MPs))

    Norm_x = np.ones_like(vx_all)
    Norm_y = np.zeros_like(vx_all)

    divider = np.sqrt(
        (vx_all * vx_all + vy_all * vy_all + vz_all * vz_all) * (Norm_x * Norm_x + Norm_y * Norm_y))
    flag_divide_by_zero = divider != 0
    cos_theta_emit = (vx_all[flag_divide_by_zero] * Norm_x[flag_divide_by_zero] + vy_all[flag_divide_by_zero] * Norm_y[flag_divide_by_zero]) / divider[flag_divide_by_zero]

    print('Done extracting emission angles.')

    return cos_theta_emit


# Test parameters
cos_theta_test = 0.7
E_impact_eV_test = 300.
n_rep = int(5e5)
# Extraction
cos_theta_emit = extract_emission_angles(n_rep, E_impact_eV_test, cos_theta_test, charge=qe, mass=me)

plt.close('all')
ms.mystyle_arial(25)

# Plotting emission angles
plt.figure(1, facecolor='white', figsize=(12, 9))
plt.hist(np.arccos(cos_theta_emit), density=True, bins=60)
plt.xlabel(r'$\theta_{emit}$')
plt.ylabel('Normalized emission angle spectrum')
plt.title('cos_theta_imp=%.2f' % cos_theta_test)

plt.figure(2, facecolor='white', figsize=(12, 9))
plt.hist(cos_theta_emit, density=True, bins=60)
plt.xlabel(r'cos($\theta_{emit}$)')
plt.title('cos_theta_imp=%.2f' % cos_theta_test)
plt.ylabel('Normalized emission angle spectrum')


plt.show()

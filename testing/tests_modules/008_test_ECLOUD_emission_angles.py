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

me = 9.10938356e-31

secondary_angle_distribution = 'cosine_3D'
sey_mod = ECL.SEY_model_ECLOUD(Emax=332., del_max=1.8848, R0=0.7, E_th=35., mufit=1.6636, secondary_angle_distribution=secondary_angle_distribution,
                               sigmafit=1.0828, switch_no_increase_energy=0, thresh_low_energy=-1)


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

    divider = np.sqrt((vx_all * vx_all + vy_all * vy_all + vz_all * vz_all) * (Norm_x * Norm_x + Norm_y * Norm_y))
    flag_divide_by_zero = divider != 0
    cos_theta_emit = (vx_all[flag_divide_by_zero] * Norm_x[flag_divide_by_zero] + vy_all[flag_divide_by_zero] * Norm_y[flag_divide_by_zero]) / divider[flag_divide_by_zero]

    print('Done extracting emission angles.')

    return cos_theta_emit, divider, flag_divide_by_zero


# Test parameters
cos_theta_test = 0.7
E_impact_eV_test = 300.
n_rep = int(5e5)
# Extraction
cos_theta_emit, divider, flag_divide_by_zero = extract_emission_angles(n_rep, E_impact_eV_test, cos_theta_test, charge=qe, mass=me)

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

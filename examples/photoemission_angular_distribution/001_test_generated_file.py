import os, sys
BIN = os.path.expanduser("../../")
sys.path.append(BIN)

import numpy as np
import matplotlib.pyplot as pl

import gen_photoemission_class as gp
import geom_impact_ellip as gie
import MP_system as mps

file_to_test = 'inv_cdf_phem_ang_distrib.mat'

# generate a chamber
chamb = gie.ellip_cham_geom_object(x_aper=2e-2, y_aper=1e-2)

# Build object used by pyecloud
ob = gp.photoemission_from_file(inv_CDF_all_photoem_file=file_to_test,
                                chamb=chamb, resc_fac=1.,
                                energy_distribution='gaussian',
                                e_pe_sigma=3.,
                                e_pe_max=10.,
                                k_pe_st=1.,
                                out_radius=10e-2,
                                photoelectron_angle_distribution='cosine_3D',
                                beamtim=None,
                                flag_continuous_emission=False)

N_mp_test = 100000
N_tests = 100

hist_list = []

for irep in range(N_tests):
    print('Test %d/%d'%(irep + 1, N_tests))
    # Build MP object
    MPe = mps.MP_system(
        N_mp_max=N_mp_test,
        nel_mp_ref_0=1.,
        fact_split=1.,
        fact_clean=1e-10,
        N_mp_regen_low=0,
        N_mp_regen=10e9,
        N_mp_after_regen=0.5e9,
        Dx_hist_reg=1, Nx_reg=10, Ny_reg=10,
        Nvx_reg=10,
        Nvy_reg=10,
        Nvz_reg=10,
        regen_hist_cut=1.,
        chamb=chamb,
        N_mp_soft_regen=10e9,
        N_mp_after_soft_regen=10e9,
        charge=-1.602176565e-19,
        mass=9.10938291e-31)

    ob.generate(MPe, lambda_t=float(N_mp_test), Dt=1. * 3e-8 * 0.1)

    theta_part = np.arctan2(MPe.y_mp[:MPe.N_mp], MPe.x_mp[:MPe.N_mp])

    hist, x_hist = np.histogram(theta_part, bins=1000, range=(-np.pi, np.pi))

    hist_list.append(hist)

combined_hist = np.mean(np.array(hist_list), axis=0)
combined_hist_norm = combined_hist / np.trapz(combined_hist, 0.5 * (x_hist[:-1] + x_hist[1:]))

pl.close('all')
pl.figure(1)
pl.plot(0.5 * (x_hist[:-1] + x_hist[1:]), combined_hist_norm)

# plot original data
input_fname = 'smooth200T_radianes.hist'
data = np.loadtxt(input_fname)
theta_data = data[:, 0]
distr_data = data[:, 1]
dist_norm = distr_data / np.trapz(distr_data, theta_data)
pl.plot(theta_data, dist_norm, 'r')

pl.show()


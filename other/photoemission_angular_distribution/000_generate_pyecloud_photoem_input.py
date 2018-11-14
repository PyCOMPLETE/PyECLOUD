import numpy as np
import matplotlib.pyplot as pl

from scipy.integrate import cumtrapz
import scipy.io as sio

# Load data
input_fname = 'smooth200T_radianes.hist'
data = np.loadtxt(input_fname)
theta_data = data[:, 0]
distr_data = data[:, 1]

# Resample data
N_theta_unif = 10000
theta_unif = np.linspace(-np.pi, np.pi, N_theta_unif)
distr_unif = np.interp(theta_unif, theta_data, distr_data)

# Build CDF
integ = cumtrapz(distr_unif, theta_unif)
cdf = np.concatenate((np.array([0.]), integ)) / np.max(integ)

# Invert
N_u_sam = 10000
u_sam = np.linspace(0, 1, N_u_sam)
angles = np.interp(u_sam, cdf, theta_unif)

# Save
pyecl_phem_filename = 'inv_cdf_phem_ang_distrib'
sio.savemat(pyecl_phem_filename, {
    'u_sam': u_sam,
    'angles': angles}, oned_as='row')


pl.close('all')
pl.figure(1)
pl.plot(theta_data, distr_data)
pl.plot(theta_unif, distr_unif, 'r')
pl.ylim(bottom=0.)

pl.figure(2)
pl.plot(theta_unif, cdf)


pl.show()


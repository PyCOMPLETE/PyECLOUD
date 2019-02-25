import sys, os
BIN=os.path.expanduser('../PyFRIENDS/')
sys.path.append(BIN)

import numpy as np
import matplotlib.pyplot as plt

import PyPIC.FiniteDifferences_ShortleyWeller_SquareGrid as PIC_FDSW
import PyPIC.geom_impact_poly as poly
import PyECLOUD.mystyle as ms

plt.close('all')
ms.mystyle_arial(fontsz=14)

x_aper = 1.5e-2
y_aper = 0.5e-2
Dh = .020e-3

x_patch_center = 1e-2 #x_aper/2.
y_patch_center = y_aper-Dh #0.
Dx_patch = 3e-3
Dy_patch = Dh/2.

Sigma_C_m2 = 1e-6

Q_tot_C_m = Sigma_C_m2*Dx_patch #1.e-12 C/mm^2

na = np.array

chamber = poly.polyg_cham_geom_object({'Vx':na([x_aper, -x_aper, -x_aper, x_aper]),
									   'Vy':na([y_aper, y_aper, -y_aper, -y_aper]),
									   'x_sem_ellip_insc':0.99*x_aper,
									   'y_sem_ellip_insc':0.99*y_aper})

pic = PIC_FDSW.FiniteDifferences_ShortleyWeller_SquareGrid(chamb = chamber, Dh = Dh, sparse_solver = 'PyKLU')


YY,XX = np.meshgrid(pic.yg, pic.xg)
rho_mat = 0*XX

mask_charge = np.logical_and(
				np.abs(XX - x_patch_center)<Dx_patch/2.,
				np.abs(YY - y_patch_center)<=Dy_patch/2.)

rho_mat[mask_charge] = 1.

rho_mat=rho_mat/(np.sum(rho_mat[:])*Dh*Dh)*Q_tot_C_m
#rho_mat[mask_charge] = np.cos(2*np.pi*XX[mask_charge]/x_aper/4.)

pic.solve(rho=rho_mat)


f1 = plt.figure(1)
f1.set_facecolor('w')
sp0 = plt.subplot(2,1,1)
plt.pcolormesh(pic.xg, pic.yg, rho_mat.T)
plt.axis('equal')
plt.grid('on')
plt.colorbar()
sp01 = plt.subplot(2,1,2, sharex=sp0)
plt.pcolormesh(pic.xg, pic.yg, pic.phi.T)
plt.axis('equal')
plt.colorbar()
plt.subplots_adjust(hspace=.3)

# plt.figure(2)
# sp1 = plt.subplot(2,1,1, sharex=sp0)
# plt.pcolormesh(pic.xg, pic.yg, pic.efx.T)
# plt.axis('equal')
# plt.colorbar()
# 
# sp2 = plt.subplot(2,1,2, sharex=sp0)
# plt.pcolormesh(pic.xg, pic.yg, pic.efy.T)
# plt.axis('equal')
# plt.colorbar()
# 
# plt.subplots_adjust(hspace=.3)

i_center = np.argmin(np.abs(pic.xg - x_patch_center))

fig3 = plt.figure()
ax1d = fig3.add_subplot(111)
ax1d.plot(pic.yg, pic.phi[i_center, :])

from scipy.constants import epsilon_0
phi_i = Sigma_C_m2*(y_aper-y_patch_center)/epsilon_0

ax1d.plot([-y_aper, y_patch_center, y_aper], [0, phi_i, 0.])

plt.show()




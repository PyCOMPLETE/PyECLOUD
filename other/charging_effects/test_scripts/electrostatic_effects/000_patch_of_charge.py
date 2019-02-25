import sys, os
BIN=os.path.expanduser('../PyFRIENDS/')
sys.path.append(BIN)

import numpy as np
import matplotlib.pyplot as plt

import PyPIC.FiniteDifferences_ShortleyWeller_SquareGrid as PIC_FDSW
import PyPIC.geom_impact_poly as poly
import PyECLOUD.mystyle as ms

plt.close('all')
ms.mystyle_arial(fontsz=16)

x_aper = 0.8e-2
y_aper = 0.5e-2
Dh = .020e-3

y_p_dist_vect = np.array([100, 50, 20, 10, 5, 2, 1 ])*Dh

x_patch_center = 0.6e-2 #x_aper/2.
Dx_patch = 3e-3
Dy_patch = Dh/2.

Sigma_C_m2 = 1e-5

Q_tot_C_m = Sigma_C_m2*Dx_patch #1.e-12 C/mm^2

na = np.array

chamber = poly.polyg_cham_geom_object({'Vx':na([x_aper, -x_aper, -x_aper, x_aper]),
    'Vy':na([y_aper, y_aper, -y_aper, -y_aper]),
    'x_sem_ellip_insc':0.99*x_aper,
    'y_sem_ellip_insc':0.99*y_aper})

pic = PIC_FDSW.FiniteDifferences_ShortleyWeller_SquareGrid(chamb = chamber, Dh = Dh, sparse_solver = 'PyKLU')


YY,XX = np.meshgrid(pic.yg, pic.xg)
max_phi_list = []
for y_patch_distance in y_p_dist_vect:
    rho_mat = 0*XX
    
    y_patch_center = y_aper - y_patch_distance
    mask_charge = np.logical_and(
    				np.abs(XX - x_patch_center)<Dx_patch/2.,
    				np.abs(YY - y_patch_center)<=Dy_patch/2.)
    
    rho_mat[mask_charge] = 1.
    
    rho_mat=rho_mat/(np.sum(rho_mat[:])*Dh*Dh)*Q_tot_C_m
    #rho_mat[mask_charge] = np.cos(2*np.pi*XX[mask_charge]/x_aper/4.)
    
    pic.solve(rho=rho_mat)
    max_phi_list.append(np.max(pic.phi))

f1 = plt.figure(1)
f1.set_facecolor('w')
sp0 = plt.subplot(1,1,1)
plt.pcolormesh(pic.xg*1e3, pic.yg*1e3, rho_mat.T)
plt.axis('equal')
plt.grid('on')
plt.colorbar()

f100 = plt.figure(100)
f100.set_facecolor('w')
sp01 = plt.subplot(1,1,1, sharex=sp0)
plt.pcolormesh(pic.xg*1e3, pic.yg*1e3, pic.phi.T)
plt.axis('equal')
cbar = plt.colorbar()
plt.subplots_adjust(bottom=.16)
cbar.set_label('Potential [V]')
sp01.set_xlim(-x_aper*1.05*1e3, x_aper*1.05*1e3) 
sp01.set_xlabel('x [mm]')
sp01.set_ylabel('y [mm]')
sp01.grid(True)

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

fig3 = plt.figure(3)
fig3.set_facecolor('w')
ax1d = fig3.add_subplot(111)
ax1d.plot(pic.yg*1e3, pic.phi[i_center, :],
    linewidth=2, label='PIC')

from scipy.constants import epsilon_0
phi_i = Sigma_C_m2*(y_aper-y_patch_center)/epsilon_0
ax1d.plot(np.array(
    [-y_aper, y_patch_center, y_aper])*1e3, [0, phi_i, 0.], 
    linewidth=2, label='Formula')
ax1d.grid(True)
ax1d.legend(loc='upper left', prop={'size':16})
ax1d.set_ylabel('Potential [V]')
ax1d.set_xlabel('y [mm]')

fig10 = plt.figure(10)
fig10.set_facecolor('w')
axll = fig10.add_subplot(1,1,1)
axll.loglog(y_p_dist_vect, max_phi_list, '.-',
        linewidth=2, markersize=10, color='b',
        label='PIC')

y_dist_theo = np.array([1e-6, 1e-2])
axll.loglog(y_dist_theo, Sigma_C_m2*y_dist_theo/epsilon_0,
        linewidth=2, color='g', label='Formula')
axll.grid(True)
axll.set_xlabel('d [m]')
axll.set_ylabel('Max. potential [V]')
axll.legend(loc='upper left', prop={'size':16})

for ff in [f1, fig3, fig10, f100]:
    ff.suptitle('Charge density: %.1e C/mm^2\nThickness:%.1e m'%(Sigma_C_m2*1e-6, y_patch_distance))
    ff.subplots_adjust(bottom=.16, top=1.-.16)

plt.show()




import sys, os
BIN=os.path.expanduser('../PyFRIENDS/')
sys.path.append(BIN)

import numpy as np
import matplotlib.pyplot as plt

import PyPIC.FiniteDifferences_ShortleyWeller_SquareGrid as PIC_FDSW
import PyPIC.geom_impact_poly as poly
import PyPIC.mystyle as ms


plt.close('all')
ms.mystyle_arial(fontsz=14)

x_aper = 10e-2
y_aper = 2e-2
Dh = .2e-3

Dx_patch = 1.
Dy_patch = 1e-3
x_patch_center = 0.
y_patch_center_vect = np.arange(y_aper-8e-3, y_aper-2e-3, 1e-3)


na = np.array

chamber = poly.polyg_cham_geom_object({'Vx':na([x_aper, -x_aper, -x_aper, x_aper]),
									   'Vy':na([y_aper, y_aper, -y_aper, -y_aper]),
									   'x_sem_ellip_insc':0.99*x_aper,
									   'y_sem_ellip_insc':0.99*y_aper})

pic = PIC_FDSW.FiniteDifferences_ShortleyWeller_SquareGrid(chamb = chamber, Dh = Dh, sparse_solver = 'PyKLU')


YY,XX = np.meshgrid(pic.yg, pic.xg)

y_probes = np.linspace(-1.1*y_aper, 1.1*y_aper, 1001)
x_probes = y_probes*0

ff = plt.figure(1)
ff.set_facecolor('w')
Ey_external = []
Ey_internal = []
for i_plot, y_patch_center in enumerate(y_patch_center_vect):
	rho_mat = 0*XX

	mask_charge = np.logical_and(
					np.abs(XX - x_patch_center)<Dx_patch/2.,
					np.abs(YY - y_patch_center)<Dy_patch/2.)

	rho_mat[mask_charge] = 1.

	# Re-normalize the charge
	rho_mat = rho_mat/np.sum(rho_mat)

	pic.solve(rho=rho_mat)

	Ex_atprb, Ey_atprb = pic.gather(x_probes, y_probes)

	color = ms.colorprog(i_plot, len(y_patch_center_vect))
	plt.plot(y_probes, Ey_atprb, color=color)


	Ey_external.append(pic.gather(na([0]), na([y_patch_center - 1.1*(Dy_patch/2)]))[1][0])
	Ey_internal.append(pic.gather(na([0]), na([y_patch_center + 1.1*(Dy_patch/2)]))[1][0])

fig1 = plt.figure(100, figsize=(8*2,6))
fig1.set_facecolor('w')
sp1 = plt.subplot(1,2,1)
sp1.plot(y_aper-y_patch_center_vect, np.abs(Ey_external), '.-')
sp2 = plt.subplot(1,2,2)
sp2.plot(y_aper-y_patch_center_vect, np.abs(Ey_internal), '.-')
for sp in [sp1,sp2]:
	sp.set_ylim(bottom=0)
	sp.set_xlabel('Distance from chamber surface [mm]')
	sp.grid('on')

sp1.set_ylabel('Electric field on the beam side')
sp2.set_ylabel('Electric field between the layer and the chamber')

fig1.subplots_adjust(bottom=.14)

plt.show()





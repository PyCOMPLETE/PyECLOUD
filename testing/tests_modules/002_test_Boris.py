import sys
from numpy import array
if '../../' not in sys.path:
    sys.path.append('../../')

import dynamics_dipole as dyndip
import dynamics_Boris_f2py as dynB

import MP_system as MPs
from geom_impact_ellip import ellip_cham_geom_object

Dt = 25e-12
N_steps = 10000
B = 0.1
N_sub_steps = 10

dynamicsd = dyndip.pusher_dipole_magnet(Dt, B)
dynamicsB = dynB.pusher_Boris(Dt, 0., B, 0.,
                              None, None, None, N_sub_steps=N_sub_steps)

chamb = ellip_cham_geom_object(.02, .02)
N_mp_max = 1000
nel_mp_ref_0 = -1
fact_split = -1
fact_clean = -1
N_mp_regen_low = -1
N_mp_regen = -1
N_mp_after_regen = -1
Dx_hist = -1
Nx_regen = -1
Ny_regen = -1
Nvx_regen = -1
Nvy_regen = -1
Nvz_regen = -1
regen_hist_cut = -1


MP_eB = MPs.MP_system(N_mp_max, nel_mp_ref_0, fact_split, fact_clean,
                      N_mp_regen_low, N_mp_regen, N_mp_after_regen,
                      Dx_hist, Nx_regen, Ny_regen, Nvx_regen, Nvy_regen, Nvz_regen, regen_hist_cut, chamb)

MP_ed = MPs.MP_system(N_mp_max, nel_mp_ref_0, fact_split, fact_clean,
                      N_mp_regen_low, N_mp_regen, N_mp_after_regen,
                      Dx_hist, Nx_regen, Ny_regen, Nvx_regen, Nvy_regen, Nvz_regen, regen_hist_cut, chamb)


N_mp = 2
Ex_n = array([0.001, 0.001])
Ey_n = array([1., 0.])

x_mpB = array([0., 0.])
y_mpB = array([0., 0.])
z_mpB = array([0., 0.])

vx_mpB = array([1., 0.])
vy_mpB = array([0., 0.])
vz_mpB = array([0., 1.])

MP_eB.x_mp = x_mpB.copy()
MP_eB.y_mp = y_mpB.copy()
MP_eB.z_mp = z_mpB.copy()

MP_eB.vx_mp = vx_mpB.copy()
MP_eB.vy_mp = vy_mpB.copy()
MP_eB.vz_mp = vz_mpB.copy()

MP_eB.N_mp = N_mp

MP_ed.x_mp = x_mpB.copy()
MP_ed.y_mp = y_mpB.copy()
MP_ed.z_mp = z_mpB.copy()

MP_ed.vx_mp = vx_mpB.copy()
MP_ed.vy_mp = vy_mpB.copy()
MP_ed.vz_mp = vz_mpB.copy()

MP_ed.N_mp = N_mp


x_lisB = []
y_lisB = []
z_lisB = []

x_lisd = []
y_lisd = []
z_lisd = []


for ii in range(N_steps):
    #print ii
    x_lisB.append(MP_eB.x_mp.copy())
    y_lisB.append(MP_eB.y_mp.copy())
    z_lisB.append(MP_eB.z_mp.copy())

    x_lisd.append(MP_ed.x_mp.copy())
    y_lisd.append(MP_ed.y_mp.copy())
    z_lisd.append(MP_ed.z_mp.copy())

    MP_eB = dynamicsB.step(MP_eB, Ex_n[0:N_mp], Ey_n[0:N_mp])

    MP_ed = dynamicsd.step(MP_ed, Ex_n[0:N_mp], Ey_n[0:N_mp])


x_lisB = array(x_lisB)
y_lisB = array(y_lisB)
z_lisB = array(z_lisB)

x_lisd = array(x_lisd)
y_lisd = array(y_lisd)
z_lisd = array(z_lisd)

import pylab as pl
pl.close('all')

for ii in range(len(x_lisB[1])):

    pl.figure(ii)
    pl.subplot(3, 1, 1)
    pl.plot(x_lisd[:, ii])
    pl.hold('on')
    pl.plot(x_lisB[:, ii], '.r')

    pl.subplot(3, 1, 2)
    pl.plot(y_lisd[:, ii])
    pl.hold('on')
    pl.plot(y_lisB[:, ii], '.r')

    pl.subplot(3, 1, 3)
    pl.plot(z_lisd[:, ii])
    pl.hold('on')
    pl.plot(z_lisB[:, ii], '.r')


pl.show()

import sys
from numpy import array

if '../../' not in sys.path:
    sys.path.append('../../')
import dynamics_Boris_f2py as dynB
import dynamics_strong_B_generalized as dyngen
import MP_system as MPs
from geom_impact_ellip import ellip_cham_geom_object
import numpy as np


#test E
Dt = 25e-13
N_steps = 100000
B = 0
N_sub_steps = 1
B_map_file = None  # 'analytic_qaudrupole_unit_grad'
fact_Bmap = 12.
N_mp = 1
Ex_n = array([1e3])
Ey_n = array([2e3])

x_mpB = array([0.0156605126754])
y_mpB = array([0.0154002281106])
z_mpB = array([1.59147686348e-05])

vx_mpB = array([335546.817425])
vy_mpB = array([-170848.411391])
vz_mpB = array([223792.460031])


#===============================================================================
# #working
# Dt=25e-13;
# N_steps=100000
# B=0
# N_sub_steps=1
# B_map_file='analytic_qaudrupole_unit_grad'
# fact_Bmap=12.*0.1
# N_mp=1
# Ex_n=array([0.,0.])
# Ey_n=array([0.,0.])
#
# x_mpB=array([0.0156605126754])
# y_mpB=array([0.0154002281106])
# z_mpB=array([1.59147686348e-05])
#
# vx_mpB=array([335546.817425])
# vy_mpB=array([-170848.411391])
# vz_mpB=array([223792.460031])
#===============================================================================


#===============================================================================
# #damped
# Dt=25e-12;
# N_steps=1000000
# B=0
# N_sub_steps=1
# B_map_file='analytic_qaudrupole_unit_grad'
# fact_Bmap=12.#*0.6
# N_mp=1
# Ex_n=array([0.,0.])
# Ey_n=array([0.,0.])
#
# x_mpB=array([0.0156605126754])
# y_mpB=array([0.0154002281106])
# z_mpB=array([1.59147686348e-05])
#
# vx_mpB=array([335546.817425])
# vy_mpB=array([-170848.411391])
# vz_mpB=array([223792.460031])
#===============================================================================

#===============================================================================
# #diverging
# Dt=25e-12;
# N_steps=100000
# B=0
# N_sub_steps=1
# B_map_file='analytic_qaudrupole_unit_grad'
# fact_Bmap=12.*2#*0.6
# N_mp=1
# Ex_n=array([0.,0.])
# Ey_n=array([0.,0.])
#
# x_mpB=array([0.0156605126754])
# y_mpB=array([0.0154002281106])
# z_mpB=array([1.59147686348e-05])
#
# vx_mpB=array([335546.817425])
# vy_mpB=array([-170848.411391])
# vz_mpB=array([223792.460031])
#===============================================================================


dynamicsB = dynB.pusher_Boris(Dt, 0., B, 0.,
                              B_map_file, fact_Bmap, None, N_sub_steps=N_sub_steps)

dynamicsGen = dyngen.pusher_strong_B_generalized(Dt, 0., B,
                                                 B_map_file, fact_Bmap, 1e-6)

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


MP_eg = MPs.MP_system(N_mp_max, nel_mp_ref_0, fact_split, fact_clean,
                      N_mp_regen_low, N_mp_regen, N_mp_after_regen,
                      Dx_hist, Nx_regen, Ny_regen, Nvx_regen, Nvy_regen, Nvz_regen, regen_hist_cut, chamb)


MP_eB.x_mp = x_mpB.copy()
MP_eB.y_mp = y_mpB.copy()
MP_eB.z_mp = z_mpB.copy()

MP_eB.vx_mp = vx_mpB.copy()
MP_eB.vy_mp = vy_mpB.copy()
MP_eB.vz_mp = vz_mpB.copy()

MP_eB.N_mp = N_mp


MP_eg.x_mp = x_mpB.copy()
MP_eg.y_mp = y_mpB.copy()
MP_eg.z_mp = z_mpB.copy()

MP_eg.vx_mp = vx_mpB.copy()
MP_eg.vy_mp = vy_mpB.copy()
MP_eg.vz_mp = vz_mpB.copy()

MP_eg.N_mp = N_mp


x_lisB = []
y_lisB = []
z_lisB = []

x_lisg = []
y_lisg = []
z_lisg = []

for ii in range(N_steps):
    if np.mod(ii, 10000) == 0:
        print(float(ii) / float(N_steps))
    #print ii
    x_lisB.append(MP_eB.x_mp.copy())
    y_lisB.append(MP_eB.y_mp.copy())
    z_lisB.append(MP_eB.z_mp.copy())

    x_lisg.append(MP_eg.x_mp.copy())
    y_lisg.append(MP_eg.y_mp.copy())
    z_lisg.append(MP_eg.z_mp.copy())

    MP_eB = dynamicsB.step(MP_eB, Ex_n[0:N_mp], Ey_n[0:N_mp])
    MP_eg = dynamicsGen.step(MP_eg, Ex_n[0:N_mp], Ey_n[0:N_mp])


x_lisB = array(x_lisB)
y_lisB = array(y_lisB)
z_lisB = array(z_lisB)

x_lisg = array(x_lisg)
y_lisg = array(y_lisg)
z_lisg = array(z_lisg)

import pylab as pl
pl.close('all')

for ii in range(len(x_lisB[1])):

    pl.figure(ii)
    sp1 = pl.subplot(3, 1, 1)
    pl.plot(x_lisB[:, ii], '.-')
    pl.plot(x_lisg[:, ii], '.r')

    pl.subplot(3, 1, 2, sharex=sp1)
    pl.plot(y_lisB[:, ii], '.-')
    pl.plot(y_lisg[:, ii], '.r')

    pl.subplot(3, 1, 3, sharex=sp1)
    pl.plot(z_lisB[:, ii], '.-')
    pl.plot(z_lisg[:, ii], '.r')
    #pl.figure(100+ii)
    #pl.plot(y_lisB[:,ii], z_lisB[:,ii],'.-')
    #pl.axis('equal')
    #pl.xlabel('y')
    #pl.ylabel('z')

    pl.figure(200 + ii)
    pl.plot(x_lisB[:, ii], y_lisB[:, ii], '.-')
    pl.plot(x_lisg[:, ii], y_lisg[:, ii], '.r')
    #pl.xlim(-np.max(np.abs(x_lisB[:,ii])), np.max(np.abs(x_lisB[:,ii])))
    #pl.ylim(-np.max(np.abs(y_lisB[:,ii])), np.max(np.abs(y_lisB[:,ii])))
    pl.axis('equal')
    pl.xlabel('x')
    pl.ylabel('y')

    #pl.figure(300+ii)
    #pl.plot(x_lisB[:,ii], z_lisB[:,ii],'.-')
    #pl.axis('equal')
    #pl.xlabel('x')
    #pl.ylabel('z')

pl.show()

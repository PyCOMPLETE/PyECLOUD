
import sys
import os
import math
import matplotlib.pyplot as plt
import numpy as np
import argparse

if '../../' not in sys.path:
    sys.path.append('../../')
import dynamics_Boris_multipole as dbu
from geom_impact_ellip import ellip_cham_geom_object
import MP_system as MPs
import mystyle as ms

from scipy.constants import m_e, e as q_e

parser = argparse.ArgumentParser()
parser.add_argument('-o', help='Path where to save the figure')
parser.add_argument('--noshow', help='Do not make a plot', action='store_true')
args = parser.parse_args()

ms.mystyle(12)
plt.close('all')

Dt = 5e-13
N_steps = 1

chamb = ellip_cham_geom_object(3., 3.)
N_mp_max = 1e4
nel_mp_ref_0 = -1
fact_split = -1
fact_clean = -1
N_mp_regen_low = -1
N_mp_regen = -1
N_mp_after_regen = -1
Dx_hist = 1.
Nx_regen = -1
Ny_regen = -1
Nvx_regen = -1
Nvy_regen = -1
Nvz_regen = -1
regen_hist_cut = -1


xx_raw = np.arange(-1, 1.01, 0.2) * 1
N_mp = int(len(xx_raw)**2)

fig2 = ms.figure('Fields at x or y = 0')
sp2 = plt.subplot(2, 2, 1)
sp2.set_title('At y = 0')
sp2.set_xlabel('x [m]')
sp2.set_ylabel(r'$B_y$ [T]')
sp2.grid(True)

sp3 = plt.subplot(2, 2, 2)
sp3.set_title('At x = 0')
sp3.set_xlabel('y [m]')
sp3.set_ylabel(r'$B_x$ [T]')
sp3.grid(True)


angles = np.pi / 4. * np.array([0, 1, 2, 3, 4], dtype=float)
for angle_ctr, angle in enumerate(angles):
    ms.figure('Angle %.1f' % angle, figsize=(16, 12))
    angle_deg = angle / np.pi * 180
    plt.suptitle('Angle = %.2f deg' % angle_deg)
    for order, title in enumerate(['Dipole', 'Quadrupole', 'Sextupole', 'Octupole']):
        B_multip = np.zeros(order + 1)
        B_skew = np.zeros(order + 1)
        B_multip[order] = 1. * np.cos(angle * (order + 1)) * math.factorial(order)
        B_skew[order] = 1. * np.sin(angle * (order + 1)) * math.factorial(order)
        if angle == 0:
            B_skew = None
        #print title, B_multip, B_skew

        pusher = dbu.pusher_Boris_multipole(Dt, B_multip=B_multip, N_sub_steps=1, B_skew=B_skew)

        MP_eB = MPs.MP_system(N_mp_max, nel_mp_ref_0, fact_split, fact_clean,
                              N_mp_regen_low, N_mp_regen, N_mp_after_regen,
                              Dx_hist, Nx_regen, Ny_regen, Nvx_regen, Nvy_regen, Nvz_regen, regen_hist_cut, chamb)

        x_mpB, y_mpB = np.meshgrid(xx_raw, xx_raw)
        x_mpB = x_mpB.flatten()
        y_mpB = y_mpB.flatten()
        z_mpB = np.zeros_like(x_mpB)

        vx_mpB = np.zeros_like(x_mpB)
        vy_mpB = np.zeros_like(x_mpB)
        vz_mpB = np.ones_like(x_mpB) * 1e5

        MP_eB.x_mp = x_mpB.copy()
        MP_eB.y_mp = y_mpB.copy()
        MP_eB.z_mp = z_mpB.copy()

        MP_eB.vx_mp = vx_mpB.copy()
        MP_eB.vy_mp = vy_mpB.copy()
        MP_eB.vz_mp = vz_mpB.copy()

        MP_eB.N_mp = N_mp

        Ex_n = np.zeros_like(x_mpB)
        Ey_n = np.zeros_like(x_mpB)

        MP_eB = pusher.step(MP_eB, Ex_n[0:N_mp], Ey_n[0:N_mp])

        del_vx = (MP_eB.vx_mp - vx_mpB)
        del_vy = (MP_eB.vy_mp - vy_mpB)

        #by = -vz_mpB * m_e * del_vx / Dt / (-q_e)
        #bx = vz_mpB * m_e * del_vy / Dt / (-q_e)

        by = -del_vx / Dt * m_e / vz_mpB / (-q_e)
        bx = del_vy / Dt * m_e / vz_mpB / (-q_e)

        sp_ctr = order + 1
        sp = plt.subplot(4, 4, sp_ctr)
        #title2 = title + ' B_multip = %s' % B_multip
        sp.set_title(title)
        sp.set_xlabel('x [m]')

        len_ = np.sqrt(bx**2 + by**2)
        if title == 'Dipole':
            outp = sp.quiver(x_mpB, y_mpB, bx / len_, by / len_, len_, pivot='middle')
            bar = plt.colorbar(outp)
            sp.set_ylabel('y [m]')
        else:
            outp = sp.quiver(x_mpB, y_mpB, bx / len_, by / len_, len_, pivot='middle')
            bar = plt.colorbar(outp)

        bar.set_ticks([])

        lim = 1.2 * np.max(xx_raw)
        sp.set_xlim(-lim, lim)
        sp.set_ylim(-lim, lim)

        mask = np.abs(y_mpB) < 1e-12
        mask2 = np.abs(x_mpB) < 1e-12

        if angle_ctr == 0:
            color = ms.colorprog(order, 4)
            xx_plot = np.linspace(-1, 1, 1e3)
            sp2.plot(x_mpB[mask], by[mask], '.', label=title, color=color)
            sp2.plot(xx_plot, xx_plot**order, color=color)
            sp3.plot(y_mpB[mask2], bx[mask2], '.', label=title, color=color)
            if order % 2 == 0:
                sp3.plot(xx_plot, np.zeros_like(xx_plot), color=color)
            else:
                sp3.plot(xx_plot, (xx_plot)**order * (-1)**((order - 1) / 2), color=color)

    sp2.set_ylim(-1.1, 1.1)
    sp3.set_ylim(-1.1, 1.1)
    #sp3.legend(loc='upper left', bbox_to_anchor=(1,1))
sp3.legend(loc='upper right', framealpha=0.5)

if args.o:
    for num in plt.get_fignums():
        fig = plt.figure(num)
        plt.suptitle('')
        if fig != fig2:
            fig.subplots_adjust(wspace=0.1)
        fig.savefig(os.path.expanduser(args.o + '_%i.png' % num))

if not args.noshow:
    plt.show()


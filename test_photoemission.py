from __future__ import division
import os
import argparse

from scipy.constants import c, e, m_e
import matplotlib.pyplot as plt
import numpy as np

import MP_system
import geom_impact_ellip
import gen_photoemission_class

import LHCMeasurementTools.mystyle as ms

plt.close('all')
ms.mystyle()


parser = argparse.ArgumentParser()
parser.add_argument('-o')
parser.add_argument('--noshow', action='store_true')
args = parser.parse_args()


def lognormal(xx, mu, sig):
    sig2 = np.sqrt(np.log(sig**2/mu**2 +1))
    mu2 = np.log(mu) - sig2**2/2
    print sig2, mu2
    fact = 1/(xx*sig2*np.sqrt(2*np.pi))*np.exp(-(np.log(xx)-mu2)**2/(2*sig2**2))
    return fact

def gauss(xx, mu, sig):
    fact = 1/(sig*np.sqrt(2*np.pi))*np.exp(-(xx-mu)**2/(2*sig**2))
    return fact

N_mp_max = int(1e5)
k_pe_st = 1/c*N_mp_max
refl_frac= 3.8e-2
nel_mp_ref = 1
alimit = 0.05

chamb = geom_impact_ellip.ellip_cham_geom_object(1.,1.)

MP_e = MP_system.MP_system(N_mp_max, nel_mp_ref, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, chamb)


phemiss = gen_photoemission_class.photoemission('unif_no_file', k_pe_st, refl_frac, 5, 7, alimit, 0.99, 0, 1.01, chamb, 0.995)

phemiss.generate(MP_e, 1, 1)


fig = ms.figure('Test Photoemission module')

xx, yy = MP_e.x_mp, MP_e.y_mp

angles = np.arctan(yy/xx)
velocities = np.array([MP_e.vx_mp, MP_e.vy_mp, MP_e.vz_mp]).T
energies = np.sum(velocities**2, axis=1)*m_e/2/e
vx, vy, vz = MP_e.vx_mp, MP_e.vy_mp, MP_e.vz_mp
#angles_v = np.pi/2 - np.arctan(vy/vx) + angles
cos_normal = -(xx*vx + yy*vy) / (np.sqrt(vx**2+vy**2+vz**2) * np.sqrt(xx**2+yy**2))
angles_v = np.arccos(cos_normal)

sin_gen = np.random.rand(len(xx))
cos_gen = 1 - sin_gen**2

sp = plt.subplot(2,2,1)
sp.grid(True)
sp.set_xlabel('x')
sp.set_xlabel('y')
sp.set_title('Positions of new MPs')
circ = plt.Circle([0,0],1, fill=False, color='black')
sp.add_artist(circ)
sp.plot(xx, yy, '.')

sp = plt.subplot(2,2,2)
sp.grid(True)
sp.set_title('Histogram of angles\nrel. to center')
sp.set_xlabel('Angle [rad]')
sp.set_yscale('log')
sp.hist(angles, bins=40, normed=True)

xx_a = np.linspace(-np.pi/2, np.pi/2, 1e5)
yy = refl_frac/np.pi + (1-refl_frac) * gauss(xx_a, 0, alimit)
sp.plot(xx_a, yy, color='g', lw=3)

sp = plt.subplot(2,2,3)
sp.grid(True)
sp.set_title('Histogram of energies')
sp.set_xlabel('Energies [eV]')
mask = energies < 30
sp.hist(energies[mask], bins=40, normed=True)
xx = np.linspace(0.5,30,1e5)
yy = lognormal(xx, 7, 5)
sp.plot(xx,yy, color='g', lw=3)

sp = plt.subplot(2,2,4)
sp.grid(True)
sp.set_title('Histogram of velocity vector\nrel. to normal')
sp.set_xlabel(r'$\theta$ [rad]')
sp.hist(angles_v, bins=40, normed=True)
sp.plot(xx_a[xx_a>0], np.cos(xx_a[xx_a>0]), color='g', lw=3)



if args.o:
    plt.suptitle('')
    plt.subplots_adjust(hspace=0.4)
    plt.savefig(os.path.expanduser(args.o))

if not args.noshow:
    plt.show()


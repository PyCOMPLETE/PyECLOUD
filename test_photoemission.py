from __future__ import division
import os
import argparse

from scipy.constants import c, e, m_e
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
import scipy.stats as stats

import MP_system
import geom_impact_ellip
import geom_impact_poly_fast_impact as gipfi
import gen_photoemission_class

import mystyle as ms

plt.close('all')
ms.mystyle()

parser = argparse.ArgumentParser()
parser.add_argument('-o', help='Define output of plots')
parser.add_argument('--noshow', action='store_true')
parser.add_argument('--energy-dist', choices=['gaussian', 'lognormal', 'lorentz'], default='gaussian')
parser.add_argument('--angle-dist', choices=['cosine_2D', 'cosine_3D'], default='cosine_3D')
args = parser.parse_args()

arr = lambda x: np.array(x, dtype=float)

def lognormal(xx, mu, sig):
    return 1/(xx*sig*np.sqrt(2*np.pi))*np.exp(-(np.log(xx)-mu)**2/(2*sig**2))

def gauss(xx, mu, sig):
    return 1/(sig*np.sqrt(2*np.pi))*np.exp(-(xx-mu)**2/(2*sig**2))

def lorentz(xx, mu, sig):
    return stats.cauchy.pdf(xx, mu, sig)

def x_angle_dist(angle_dist,x):
    if angle_dist == 'cosine_2D':
        return np.cos(x)
    elif angle_dist == 'cosine_3D':
        return np.sin(2*x)


# Config

# Sample size of generated photoelectrons
k_pe_st = 1
refl_frac= 20e-2
nel_mp_ref = 1
alimit = 0.05

if args.energy_dist == 'lorentz':
    mu, sig = 0.64, 3.7
elif args.energy_dist == 'lognormal':
    sig = np.sqrt(np.log(5**2/7**2 +1))
    mu = np.log(7) - sig**2/2
else:
    mu, sig = 7, 5

my_chamb_dict = {
    'Vx': arr([-1,-1,1,1,1.1,1.1,1,1]),
    'Vy': arr([1,-1,-1,-0.3,-0.5,0.5,0.3,1]),
}
my_chamb_dict['x_sem_ellip_insc'] = 0.98*my_chamb_dict['Vx'].max()
my_chamb_dict['y_sem_ellip_insc'] = 0.98*my_chamb_dict['Vy'].max()

my_chamb_dict['phem_cdf'] = np.round(np.cumsum(arr([0.1,0.1, 0.1, 0.1, 0.3,0.1,0.1,0.1])), 3)

chamb_seg = gipfi.polyg_cham_geom_object(my_chamb_dict, False, flag_assume_convex=False, distance_new_phem=1e-10, flag_counter_clockwise_chamb=False)

phem_seg = gen_photoemission_class.photoemission_per_segment(chamb_seg, args.energy_dist, sig, mu, k_pe_st, args.angle_dist)



chamb_elip = geom_impact_ellip.ellip_cham_geom_object(1.,1.)
phem_elip = gen_photoemission_class.photoemission('unif_no_file', k_pe_st, refl_frac, sig, mu, alimit, 0.99, 0, 1.01, chamb_elip, 0.995, args.energy_dist, args.angle_dist, None, False)


for phem, chamb, N_mp_max in [
    (phem_elip, chamb_elip, int(1e6)),
    (phem_seg, chamb_seg, int(1e3)),
]:

    MP_e = MP_system.MP_system(N_mp_max, nel_mp_ref, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, chamb)
    phem.generate(MP_e, 1/c*N_mp_max, 1)

    xx, yy = MP_e.x_mp, MP_e.y_mp

    angles = np.arctan(yy/xx)
    velocities = np.array([MP_e.vx_mp, MP_e.vy_mp, MP_e.vz_mp]).T
    energies = np.sum(velocities**2, axis=1)*m_e/2/e
    vx, vy, vz = MP_e.vx_mp, MP_e.vy_mp, MP_e.vz_mp
    cos_normal = -(xx*vx + yy*vy) / (np.sqrt(vx**2+vy**2+vz**2) * np.sqrt(xx**2+yy**2))
    angles_v = np.arccos(cos_normal)

    sin_gen = np.random.rand(len(xx))
    cos_gen = 1 - sin_gen**2

    fig = ms.figure('Test Photoemission module %s' % phem.__class__)

    sp = plt.subplot(2,2,1)
    sp.grid(True)
    sp.set_xlabel('x')
    sp.set_ylabel('y')
    sp.set_title('Positions of %.1e new MPs' % len(xx))
    if phem is phem_elip:
        circ = plt.Circle([0,0],1, fill=False, color='black')
        sp.add_artist(circ)
    else:
        sp.plot(chamb.Vx, chamb.Vy)
        sp.set_xlim(1.1*min(chamb.Vx), 1.1*max(chamb.Vx))
        sp.set_ylim(1.1*min(chamb.Vy), 1.1*max(chamb.Vy))
    sp.plot(xx, yy, '.')

    if phem is phem_elip:
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
    if args.energy_dist == 'lognormal':
        yy = lognormal(xx, mu, sig)
    elif args.energy_dist == 'gaussian':
        factor = integrate.quad(lambda x: gauss(x, mu, sig), -np.inf, np.inf)[0] / integrate.quad(lambda x: gauss(x, mu, sig), 0, np.inf)[0]
        yy = gauss(xx, mu, sig) * factor
    elif args.energy_dist == 'lorentz':
        yy = lorentz(xx, mu, sig) / (1-stats.cauchy.cdf(0, loc=mu, scale=sig))

    sp.plot(xx,yy, color='g', lw=3)

    if phem is phem_elip:
        sp = plt.subplot(2,2,4)
        sp.grid(True)
        sp.set_title('Histogram of velocity vector\nrel. to normal')
        sp.set_xlabel(r'$\theta$ [rad]')
        sp.hist(angles_v, bins=40, normed=True)
        sp.plot(xx_a[xx_a>0], x_angle_dist(args.angle_dist, xx_a[xx_a>0]), color='g', lw=3)

        xx = my_chamb_dict['Vx'].squeeze()
        yy = my_chamb_dict['Vy'].squeeze()
        xx = list(xx) + [xx[0]]
        yy = list(yy) + [yy[0]]


if args.o:
    for num in plt.get_fignums():
        fig = plt.figure(num)
        plt.suptitle('')
        plt.subplots_adjust(hspace=0.4)
        plt.savefig(os.path.expanduser(args.o+'_%i.png' % num))

if not args.noshow:
    plt.show()


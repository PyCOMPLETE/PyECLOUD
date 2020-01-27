
import sys
import os
import argparse

from scipy.constants import c, e, m_e
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
import scipy.stats as stats

if '../../' not in sys.path:
    sys.path.append('../../')
import MP_system
import geom_impact_ellip
import geom_impact_rect_fast_impact as girfi
import geom_impact_poly_fast_impact as gipfi
import gen_photoemission_class

import mystyle as ms

plt.close('all')
ms.mystyle()

parser = argparse.ArgumentParser()
parser.add_argument('-o', help='Define output of plots')
parser.add_argument('--noshow', action='store_true')
parser.add_argument('--energy-dist', choices=['gaussian', 'lognormal', 'lorentz', 'rect', ], default='gaussian')
parser.add_argument('--angle-dist', choices=['cosine_2D', 'cosine_3D'], default='cosine_3D')
args = parser.parse_args()

arr = lambda x: np.array(x, dtype=float)


def lognormal(xx, mu, sig):
    return 1 / (xx * sig * np.sqrt(2 * np.pi)) * np.exp(-(np.log(xx) - mu)**2 / (2 * sig**2))


def gauss(xx, mu, sig):
    return 1 / (sig * np.sqrt(2 * np.pi)) * np.exp(-(xx - mu)**2 / (2 * sig**2))


def lorentz(xx, mu, sig):
    return stats.cauchy.pdf(xx, mu, sig)


def x_angle_dist(angle_dist, x):
    if angle_dist == 'cosine_2D':
        return np.cos(x)
    elif angle_dist == 'cosine_3D':
        return np.sin(2 * x)


# Config

# Sample size of generated photoelectrons
k_pe_st = 1
refl_frac = 20e-2


def init_mp(N_mp_max, chamb):
    nel_mp_ref = 1
    return MP_system.MP_system(N_mp_max, nel_mp_ref, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, chamb)

# Reasonable mu, sigma values for different energy distributions
if args.energy_dist == 'lorentz':
    mu, sig = 0.64, 3.7
elif args.energy_dist == 'lognormal':
    sig = np.sqrt(np.log(5**2 / 7**2 + 1))
    mu = np.log(7) - sig**2 / 2
else:
    mu, sig = 7, 5

# Have local variables for each module, to not mix them up

## Photoemission model '1' (traditional)
# This test is for a circular chamber.
# Here it is easy to test the angular emission distribution of the photoelectrons,
# since the normal vector of the surface is the particle position.
# It is also easy to show the theoretical function of the spatial distribution together
# with a histogram of the actual positions using the module.


def test_model_1():

    alimit = 0.05
    # Circular chamber with radius 1
    radius = 1.
    chamb_elip = geom_impact_ellip.ellip_cham_geom_object(radius, radius)

    # Photoemission object
    phem_elip = gen_photoemission_class.photoemission(
        'unif_no_file', k_pe_st, refl_frac, sig, mu, alimit, 0.99, 0, 1.01, chamb_elip, 0.995,
        args.energy_dist, args.angle_dist, None, False)

    # 1e6 number of MPs
    N_mp_max = int(1e6)

    MP_e = init_mp(N_mp_max, chamb_elip)
    phem_elip.generate(MP_e, 1 / c * N_mp_max, 1)

    xx, yy = MP_e.x_mp, MP_e.y_mp

    angles = np.arctan(yy / xx)
    velocities = np.array([MP_e.vx_mp, MP_e.vy_mp, MP_e.vz_mp]).T
    energies = np.sum(velocities**2, axis=1) * m_e / 2 / e

    sp = plt.subplot(2, 2, 1)
    sp.grid(True)
    sp.set_xlabel('X dimension')
    sp.set_ylabel('Y dimension')
    sp.set_title('Positions of %.1e new MPs' % len(xx))
    circ = plt.Circle([0, 0], 1, fill=False, color='black')
    sp.add_artist(circ)
    sp.plot(xx, yy, '.')

    sp = plt.subplot(2, 2, 2)
    sp.grid(True)
    sp.set_title('Histogram of angles\nrel. to center')
    sp.set_xlabel('Angles [rad]')
    sp.set_yscale('log')
    sp.hist(angles, bins=40, normed=True)

    xx_a = np.linspace(-np.pi / 2, np.pi / 2, 1e5)
    yy = refl_frac / np.pi + (1 - refl_frac) * gauss(xx_a, 0, alimit)
    sp.plot(xx_a, yy, color='g', lw=3)

    sp = plt.subplot(2, 2, 3)
    sp.grid(True)
    sp.set_title('Histogram of energies')
    sp.set_xlabel('Energies [eV]')
    mask = energies < 30
    sp.hist(energies[mask], bins=40, normed=True)
    xx = np.linspace(0.5, 30, 1e5)
    if args.energy_dist == 'lognormal':
        yy = lognormal(xx, mu, sig)
    elif args.energy_dist == 'gaussian':
        factor = integrate.quad(lambda x: gauss(x, mu, sig), -np.inf, np.inf)[0] / integrate.quad(lambda x: gauss(x, mu, sig), 0, np.inf)[0]
        yy = gauss(xx, mu, sig) * factor
    elif args.energy_dist == 'lorentz':
        yy = lorentz(xx, mu, sig) / (1 - stats.cauchy.cdf(0, loc=mu, scale=sig))
    elif args.energy_dist == 'rect':
        yy = np.zeros_like(xx)
        mask = np.logical_and(xx > mu - sig / 2, mu + sig / 2 > xx)
        yy[mask] = 1 / sig

    sp.plot(xx, yy, color='g', lw=3)

    xx, yy = MP_e.x_mp, MP_e.y_mp

    angles = np.arctan(yy / xx)
    velocities = np.array([MP_e.vx_mp, MP_e.vy_mp, MP_e.vz_mp]).T
    energies = np.sum(velocities**2, axis=1) * m_e / 2 / e
    vx, vy, vz = MP_e.vx_mp, MP_e.vy_mp, MP_e.vz_mp
    cos_normal = -(xx * vx + yy * vy) / (np.sqrt(vx**2 + vy**2 + vz**2) * np.sqrt(xx**2 + yy**2))
    angles_v = np.arccos(cos_normal)

    sp = plt.subplot(2, 2, 4)
    sp.grid(True)
    sp.set_title('Histogram of velocity vector\nrel. to normal')
    sp.set_xlabel(r'$\theta$ [rad]')
    sp.hist(angles_v, bins=40, normed=True)
    sp.plot(xx_a[xx_a > 0], x_angle_dist(args.angle_dist, xx_a[xx_a > 0]), color='g', lw=3)

# Photoemission model 'from_file' (2)


def test_model_2():
    # Generate a sinusoidal-squared distribution on a square chamber.
    # This also tests the rectangular chamber module
    chamber_width = 1.
    chamb_rect = girfi.rect_cham_geom_object(chamber_width, chamber_width, False)
    N_mp_max = int(1e6)

    n_dist = int(1e3)
    angles = np.linspace(-np.pi, np.pi, n_dist)
    dist = np.cumsum(np.sin(angles)**2)
    dist /= dist.max()
    dist_dict = {
        'angles': angles,
        'u_sam': dist,
    }

    phem_rect = gen_photoemission_class.photoemission_from_file(
        dist_dict, chamb_rect, 0.99, args.energy_dist, sig, mu, k_pe_st, 1.01 * np.sqrt(2), args.angle_dist, None, None)

    MP_e = init_mp(N_mp_max, chamb_rect)
    phem_rect.generate(MP_e, 1 / c * N_mp_max, 1)

    sp = plt.subplot(2, 2, 1)
    sp.grid(True)
    sp.set_xlabel('X dimension')
    sp.set_ylabel('Y dimension')
    sp.set_title('Positions of %.1e new MPs' % len(MP_e.x_mp))
    sp.plot(chamb_rect.Vx, chamb_rect.Vy)
    sp.plot(MP_e.x_mp, MP_e.y_mp, ls='None', marker='.')

    angles_generated = np.arctan2(MP_e.y_mp, MP_e.x_mp)
    hist, bins = np.histogram(angles_generated, n_dist)
    factor_hist = np.trapz(hist, angles)

    sp = plt.subplot(2, 2, 2)
    sp.grid(True)
    sp.set_xlabel('Angles [rad]')
    sp.set_ylabel('Normalized # generated')
    sp.set_title('Sin$^2$ distribution - from file')
    sp.plot(angles, np.sin(angles)**2 / np.pi, color='g', lw=3)
    sp.step(angles, hist / factor_hist, color='b')

    chamb_segment_x = np.concatenate([
        (np.linspace(-1, 1, 100))[:-1],
        (np.ones(100, float))[:-1],
        (np.linspace(1, -1, 100))[:-1],
        (np.ones(100, float) * -1)[:-1],
    ])

    chamb_segment_y = np.concatenate([
        (np.ones(100, float) * -1)[:-1],
        (np.linspace(-1, 1, 100))[:-1],
        (np.ones(100, float))[:-1],
        (np.linspace(1, -1, 100))[:-1],
    ])

    phem_pdf = []
    for seg_ctr, (seg_x, seg_y) in enumerate(zip(chamb_segment_x, chamb_segment_y)):
        angle_0 = np.arctan2(seg_y, seg_x) % np.pi

        seg_ctr_1 = (seg_ctr + 1) % len(chamb_segment_x)
        angle_1 = np.arctan2(chamb_segment_y[seg_ctr_1], chamb_segment_x[seg_ctr_1]) % np.pi

        if angle_0 > angle_1:
            angle_0 -= np.pi
        phem_pdf.append(integrate.quad(lambda x: np.sin(x)**2 / np.pi, angle_0, angle_1)[0])

    phem_cdf = np.cumsum(phem_pdf, dtype=float)
    phem_cdf /= phem_cdf[-1]

    phem_chamb_dict = {
        'Vx': chamb_segment_x,
        'Vy': chamb_segment_y,
        'x_sem_ellip_insc': np.sqrt(2) * 1.01,
        'y_sem_ellip_insc': np.sqrt(2) * 1.01,
        'phem_cdf': phem_cdf,
    }

    chamb_phem = gipfi.polyg_cham_photoemission(phem_chamb_dict)
    if chamb_phem.vertexes_are_subset(chamb_rect):
        print('Vertexes are compatible -> Good!')
    else:
        raise ValueError('Vertexes of the two chambers are not compatible!')

    phem_segment = gen_photoemission_class.photoemission_per_segment(
        chamb_phem, args.energy_dist, sig, mu, k_pe_st, args.angle_dist, None, None)

    sp = plt.subplot(2, 2, 3)
    sp.grid(True)
    sp.set_xlabel('X dimension')
    sp.set_ylabel('Y dimension')
    sp.set_title('Chamber vertexes')
    sp.set_ylim(-1.1, 1.1)
    sp.set_xlim(-1.1, 1.1)
    sp.plot(chamb_rect.Vx, chamb_rect.Vy, lw=3, marker='o', label='Rect')
    sp.plot(chamb_phem.Vx, chamb_phem.Vy, ls='None', marker='.', label='Photoem.')
    sp.legend(title='Chamber', loc='best')

    MP_e = init_mp(N_mp_max, chamb_rect)
    phem_segment.generate(MP_e, 1 / c * N_mp_max, 1)

    sp = plt.subplot(2, 2, 4)
    sp.grid(True)
    sp.set_xlabel('Angles [rad]')
    sp.set_ylabel('Normalized # generated')
    sp.set_title('Sin$^2$ distribution - per segment')

    angles_generated = np.arctan2(MP_e.y_mp, MP_e.x_mp)
    hist, bins = np.histogram(angles_generated, n_dist)
    factor_hist = np.trapz(hist, angles)

    sp.plot(angles, np.sin(angles)**2 / np.pi, color='g', lw=3)
    sp.step(angles, hist / factor_hist, color='b')

## Photoemission model 'per_segment' (3)


def test_model_3():
    # Non-convex chamber
    my_chamb_dict = {
        'Vx': arr([-1, -1, 1, 1, 1.1, 1.1, 1, 1]),
        'Vy': arr([1, -1, -1, -0.3, -0.5, 0.5, 0.3, 1]),
    }
    my_chamb_dict['x_sem_ellip_insc'] = 0.98 * my_chamb_dict['Vx'].max()
    my_chamb_dict['y_sem_ellip_insc'] = 0.98 * my_chamb_dict['Vy'].max()

    my_chamb_dict['phem_cdf'] = np.round(np.cumsum(arr([0.1, 0.1, 0.1, 0.1, 0.3, 0.1, 0.1, 0.1])), 3)

    chamb_seg = gipfi.polyg_cham_geom_object(my_chamb_dict, False, flag_assume_convex=False, flag_verbose_stdout=True)
    chamb_photo_seg = gipfi.polyg_cham_photoemission(my_chamb_dict)

    phem_seg = gen_photoemission_class.photoemission_per_segment(
        chamb_photo_seg, args.energy_dist, sig, mu, k_pe_st, args.angle_dist)
    N_mp_max = int(1e3)

    MP_e = init_mp(N_mp_max, chamb_seg)
    phem_seg.generate(MP_e, 1 / c * N_mp_max, 1)

    sp = plt.subplot(2, 2, 1)
    sp.grid(True)
    sp.set_xlabel('x')
    sp.set_ylabel('y')
    sp.set_title('Positions of %.1e new MPs' % N_mp_max)
    sp.plot(chamb_seg.Vx, chamb_seg.Vy)
    sp.set_xlim(1.1 * np.min(chamb_seg.Vx), 1.1 * np.max(chamb_seg.Vx))
    sp.set_ylim(1.1 * np.min(chamb_seg.Vy), 1.1 * np.max(chamb_seg.Vy))
    sp.plot(MP_e.x_mp, MP_e.y_mp, '.')

# Call tests
ms.figure("Test Photoemission module 'photemission'")
test_model_1()
ms.figure("Test Photoemission module 'photemission_per_segment'")
test_model_3()
ms.figure("Test Photoemission module 'photemission_from_file' and 'per_segment'")
test_model_2()

if args.o:
    for num in plt.get_fignums():
        fig = plt.figure(num)
        plt.suptitle('')
        plt.subplots_adjust(hspace=0.4)
        plt.savefig(os.path.expanduser(args.o + '_%i.png' % num))

if not args.noshow:
    plt.show()

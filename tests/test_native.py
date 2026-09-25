"""Checks for the installed Fortran and Cython extensions."""

from importlib import import_module

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
import pytest
from scipy.special import wofz


@pytest.mark.parametrize('name', [
    'rhocompute', 'int_field_for', 'hist_for', 'seg_impact', 'errffor',
    'boris_step', 'vectsum', 'boris_cython', 'geom_impact_poly_cython',
    'buildup_simulation', 'PyEC4PyHT',
])
def test_imports(name):
    import_module('PyECLOUD.' + name)


def test_charge_deposition_and_field_interpolation():
    from PyECLOUD import rhocompute, int_field_for

    x, y = np.array([0.25, 1.5]), np.array([0.5, 1.25])
    weights = np.array([2., 3.])
    rho = rhocompute.compute_sc_rho(x, y, weights, 0., 0., 1., 4, 4)
    assert_allclose(rho.sum(), weights.sum())
    assert_allclose(rho[0, 0], 2 * 0.75 * 0.5)
    xx, yy = np.meshgrid(np.arange(4.), np.arange(4.), indexing='ij')
    ex, ey = int_field_for.int_field(x, y, 0., 0., 1., 1., 2*xx + yy, xx - yy)
    assert_allclose(ex, 2*x + y)
    assert_allclose(ey, x - y)


def test_histograms_and_sum():
    from PyECLOUD import hist_for, seg_impact, vectsum

    weights = np.array([2., 4., 3.])
    hist = np.zeros(4)
    hist_for.compute_hist(np.array([0.25, 1.5, 2.]), weights, 0., 1., hist)
    assert_allclose(hist, [1.5, 2.5, 5., 0.])
    assert_allclose(vectsum.vectsum(weights), 9.)
    impacts = np.zeros(3)
    seg_impact.update_seg_impact(np.array([0, 2, 2, -1, 3], dtype=np.int32),
                                  np.array([1., 2., 3., 10., 10.]), impacts)
    assert_allclose(impacts, [1., 0., 5.])


def test_complex_error_function():
    from PyECLOUD.errffor import errf

    for z in (0j, 0.3 + 0.7j, 2 + 1j):
        real, imag = errf(z.real, z.imag)
        assert_allclose(real + 1j*imag, wofz(z), rtol=2e-6, atol=1e-10)


@pytest.mark.parametrize('backend', ['fortran', 'cython'])
def test_boris_constant_electric_field(backend):
    from PyECLOUD.boris_step import boris_step
    from PyECLOUD.boris_cython import boris_step_multipole

    state = [np.zeros(2) for _ in range(6)]
    field = [np.array([1., 2.]), np.array([-2., 3.]), np.zeros(2)]
    magnetic = [np.zeros(2) for _ in range(3)]
    dt, charge, mass = 0.1, -1., 2.
    if backend == 'fortran':
        boris_step(dt, *state, *field, *magnetic, mass, charge)
    else:
        boris_step_multipole(1, dt, np.zeros(1), np.zeros(1),
                            *state, *field[:2], *magnetic, True, charge, mass)
    expected_velocity = charge / mass * dt * np.array(field)
    assert_allclose(state[3:], expected_velocity)
    assert_allclose(state[:3], dt * expected_velocity)


def test_boris_backends_agree_in_magnetic_field():
    from PyECLOUD.boris_step import boris_step
    from PyECLOUD.boris_cython import boris_step_multipole

    state_f = [np.zeros(2) for _ in range(3)] + [np.array([1., 2.]) for _ in range(3)]
    state_c = [arr.copy() for arr in state_f]
    initial_speed_squared = np.sum(np.array(state_f[3:])**2, axis=0)
    zero = np.zeros(2)
    magnetic = [np.full(2, 0.2), np.full(2, 0.3), np.full(2, 0.4)]
    for _ in range(20):
        boris_step(0.1, *state_f, zero, zero, zero, *magnetic, 2., -1.)
    boris_step_multipole(20, 0.1, np.zeros(1), np.zeros(1),
                        *state_c, zero, zero, *magnetic, True, -1., 2.)
    assert_allclose(state_c, state_f, rtol=1e-13, atol=1e-14)
    assert_allclose(np.sum(np.array(state_c[3:])**2, axis=0), initial_speed_squared)


def test_polygon_classification_and_impact():
    from PyECLOUD.geom_impact_rect_fast_impact import rect_cham_geom_object

    chamber = rect_cham_geom_object(1., 0.5, flag_non_unif_sey=False)
    assert_array_equal(chamber.is_outside(np.array([0., 1.1, 0.]),
                                        np.array([0., 0., 0.6])), [False, True, True])
    zero = np.zeros(1)
    x, y, z, nx, ny, _ = chamber.impact_point_and_normal(
        zero, zero, zero, np.array([2.]), zero, zero)
    assert_allclose(x, [0.99])
    assert_allclose(y, 0.)
    assert_allclose(z, 0.)
    assert_allclose(np.abs(nx), 1.)
    assert_allclose(ny, 0.)

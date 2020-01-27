
import numpy as np
import scipy.io as sio


def get_linear_extrapolation_parameters(x, y, range_extrapolate_right):
    # Find fit parameters for energies that extend the measured ranges
    mask_xx_fit = x > (x.max() - range_extrapolate_right)
    xx_fit = x[mask_xx_fit]
    if len(xx_fit) < 2:
        raise ValueError('Range for SEY fit is too small! You may have to increase the range_extrapolate_right parameter!')
    yy_fit = y[mask_xx_fit]
    extrapolate_grad, extrapolate_const = np.polyfit(xx_fit, yy_fit, 1)

    return extrapolate_grad, extrapolate_const


def resample_sey_data(energy_eV_samples, sey_true_samples, sey_elast_samples, uniform_dE, range_extrapolate_right):

    energy_eV_0 = np.array(energy_eV_samples, dtype=float)
    sey_true_0 = np.array(sey_true_samples, dtype=float)
    sey_elast_0 = np.array(sey_elast_samples, dtype=float)

    # Build equally spaced arrays that are used by the interp function
    energy_eV = np.arange(energy_eV_0.min(), energy_eV_0.max() + uniform_dE * .5, uniform_dE)

    sey_true = np.interp(energy_eV, energy_eV_0, sey_true_0)
    extrapolate_grad_true, extrapolate_const_true = get_linear_extrapolation_parameters(energy_eV, sey_true, range_extrapolate_right)

    sey_elast = np.interp(energy_eV, energy_eV_0, sey_elast_0)
    extrapolate_grad_elast, extrapolate_const_elast = get_linear_extrapolation_parameters(energy_eV, sey_elast, range_extrapolate_right)

    resampled = {
        'energy_eV': energy_eV,
        'sey_true': sey_true,
        'sey_elast': sey_elast,
        'extrapolate_grad_true': extrapolate_grad_true,
        'extrapolate_const_true': extrapolate_const_true,
        'extrapolate_grad_elast': extrapolate_grad_elast,
        'extrapolate_const_elast': extrapolate_const_elast,
    }

    return resampled



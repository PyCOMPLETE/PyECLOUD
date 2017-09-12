from __future__ import division, print_function
import os

import numpy as np
import scipy.io as sio

class SEY_model_from_file(object):

    def __init__(self, sey_file, flag_factor_costheta):
        """
        - sey file is the path to a mat file of the correct format, either an absolute path or in the sey_files folder.
        - if flag_factor_costheta is True, the SEY is increased depending on the angle with which the electrons are hitting.
            Set to None or False to disable.
        """

        # Search for file with given name, either as an absolute path or in the dedicated folder
        # Absolute path has precedence.
        candidate_files = [
            os.path.expanduser(sey_file),
            os.path.abspath(os.path.dirname(__file__)) + '/sey_files/' + sey_file,
        ]
        existing_files = filter(os.path.isfile, candidate_files)
        if not existing_files:
            raise ValueError('SEY file %s is not found' % sey_file)
        sey_file_real = existing_files[0]
        print('Secondary emission from file %s' % sey_file_real)

        sey_properties = sio.loadmat(sey_file)
        energy_eV = sey_properties['energy_eV'].squeeze()
        sey_parameter = sey_properties['sey_parameter'].squeeze()
        extrapolate_grad = float(sey_properties['extrapolate_grad'].squeeze())
        extrapolate_const = float(sey_properties['extrapolate_const'].squeeze())

        diff_e = np.round(np.diff(energy_eV), 3)
        delta_e = diff_e[0]
        if np.any(diff_e != delta_e):
            raise ValueError('Energy in file %s is not equally spaced.' % sey_file_real)

        # sey_diff is needed by the interp function
        # A 0 is appended because this last element is never needed but the array must have the correct shape
        self.sey_diff = np.append(np.diff(sey_parameter), 0.)

        # This merely populates the object namespace
        self.energy_eV              = energy_eV
        self.sey_parameter          = sey_parameter
        self.sey_file               = sey_file
        self.flag_factor_costheta   = flag_factor_costheta
        self.energy_eV_min          = energy_eV.min()
        self.energy_eV_max          = energy_eV.max()
        self.delta_e                = delta_e
        self.extrapolate_grad       = extrapolate_grad
        self.extrapolate_const      = extrapolate_const

    def SEY_process(self,nel_impact,E_impact_eV, costheta_impact, i_impact):

        delta = np.zeros_like(E_impact_eV, dtype=float)

        mask_fit = (E_impact_eV > self.energy_eV_max)
        mask_regular = ~mask_fit

        delta[mask_regular] = self.interp(E_impact_eV[mask_regular])
        delta[mask_fit] = self.extrapolate_const + self.extrapolate_grad * E_impact_eV[mask_fit]

        if self.flag_factor_costheta:
            factor_costheta = np.exp(0.5*(1.-costheta_impact))
        else:
            factor_costheta = 1.

        # The concept of reflected electrons is not honored
        ref_frac = np.zeros_like(E_impact_eV, dtype=bool)

        return nel_impact*delta*factor_costheta, ref_frac, ~ref_frac

    def interp(self, energy_eV):
        """
        Linear interpolation of the energy - SEY curve.
        """
        index_float = (energy_eV - self.energy_eV_min)/self.delta_e
        index_remainder, index_int = np.modf(index_float)
        index_int = index_int.astype(int)
        return self.sey_parameter[index_int] + index_remainder*self.sey_diff[index_int]

    def interp_regular(self, energy_eV):
        #This fails if the input is not in ascending order.
        #return np.interp(energy_eV, self.energy_eV, self.sey_parameter)
        raise ValueError('Warning! Do not use interp_regular!')


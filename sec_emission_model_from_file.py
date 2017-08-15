from __future__ import division
import os
import numpy as np

class SEY_model_from_file(object):

    def __init__(self, sey_file, range_extrapolate_right, delta_e, flag_factor_costheta, factor_sey):
        """
        - sey file is the path to a file of the correct format, either an absolute path or in the sey_files/SEY-LE_SEY folder.
            See the example files in sey_files/SEY-LE_SEY for the format that should be used for input files.
        - range_extrapolate_right states the range in electron volts which is used to create a linear fit in order to extrapolate high energies.
        - delta_e is the resolution used for the interpolation, for example .1 eV.
        - if flag_factor_costheta is True, the SEY is increased depending on the angle with which the electrons are hitting.
        - factor_sey is a factor applied to the SEY from the file, for example to emulate scrubbing effects

        """

        # Search for file with given name, either as an absolute path or in the dedicated folder
        candidate_files = [
            os.path.expanduser(sey_file),
            os.path.abspath(os.path.dirname(__file__)) + '/sey_files/SEY-LE_SEY/' + sey_file,
        ]
        existing_files = filter(os.path.isfile, candidate_files)
        if not existing_files:
            raise ValueError('SEY file %s is not found' % sey_file)
        elif len(existing_files) != 1:
            raise ValueError('More than one sey file found for given name %s: %r' % (sey_file, existing_files))
        sey_file_real = existing_files[0]

        # Parse file
        energy_eV_list = []
        sey_parameter_list = []
        with open(sey_file_real) as f:
            for ctr, line in enumerate(f):
                split = line.split()
                if line.startswith('#'):
                    continue
                elif len(split) == 2:
                    energy_eV_list.append(float(split[0]))
                    sey_parameter_list.append(float(split[1]))
        energy_eV_0 = np.array(energy_eV_list, dtype=float)
        sey_parameter_0 = np.array(sey_parameter_list, dtype=float) * factor_sey

        # Build equally spaced arrays that are used by the interp function
        energy_eV = np.arange(energy_eV_0.min(), energy_eV_0.max()+delta_e*.5, delta_e)
        sey_parameter = np.interp(energy_eV, energy_eV_0, sey_parameter_0)

        # Find fit parameters for energies that extend the measured ranges
        mask_xx_fit = energy_eV_0 > (energy_eV.max() - range_extrapolate_right)
        xx_fit = energy_eV_0[mask_xx_fit]
        if len(xx_fit) < 2:
            raise ValueError('Range for SEY fit is too small! You may have to increase the range_extrapolate_right parameter!')
        yy_fit = sey_parameter_0[mask_xx_fit]
        self.extrapolate_grad, self.extrapolate_const = np.polyfit(xx_fit, yy_fit, 1)

        # sey_diff is needed by the interp function
        # A 0 is appended because this last element is never needed but the array must have the correct shape
        self.sey_diff = np.append(np.diff(sey_parameter), 0.)

        # This merely populates the object namespace
        self.energy_eV              = energy_eV
        self.sey_parameter          = sey_parameter
        self.energy_eV_0            = energy_eV_0
        self.sey_parameter_0        = sey_parameter_0
        self.sey_file               = sey_file
        self.flag_factor_costheta   = flag_factor_costheta
        self.energy_eV_min          = energy_eV.min()
        self.energy_eV_max          = energy_eV.max()
        self.delta_e                = delta_e
        self.factor_sey             = factor_sey

        print('Secondary emission from file %s' % sey_file_real)

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


from __future__ import division
import os
import numpy as np

class SEY_from_file(object):

    def __init__(self, sey_file, R0=None, default_work_function=4.25, range_extrapolate_right=None, delta_e=0.1):
        """
        range_extrapolate_right states the range in electron volts which is used to create a linear fit in order to extrapolate high energies.
        """

        work_function = default_work_function
        energy_eV_list = []
        sey_parameter_list = []
        with open(os.path.expanduser(sey_file)) as f:
            for ctr, line in enumerate(f):
                split = line.split()
                if 'Work function' in line:
                    work_function = float(split[-1])
                elif len(split) == 2:
                    energy_eV_list.append(split[0])
                    sey_parameter_list.append(split[1])

        energy_eV_0 = np.array(energy_eV_list, dtype=float) + work_function
        sey_parameter_0 = np.array(sey_parameter_list, dtype=float)

        # Build equally spaced arrays that are used by the interp function
        energy_eV = np.arange(energy_eV_0.min(), energy_eV_0.max()+delta_e*.1, delta_e)
        sey_parameter = np.interp(energy_eV, energy_eV_0, sey_parameter_0)

        if range_extrapolate_right is not None:
            mask_xx_fit = energy_eV > energy_eV.max() - range_extrapolate_right
            xx_fit = energy_eV[mask_xx_fit]
            if len(xx_fit) < 2:
                raise ValueError('Range for fit is too small!')
            yy_fit = sey_parameter[mask_xx_fit]
            self.extrapolate_grad, self.extrapolate_const = np.polyfit(xx_fit, yy_fit, 1)
        else:
            self.extrapolate_grad, self.extrapolate_const = None, None

        self.energy_eV = energy_eV
        self.sey_parameter = sey_parameter
        self.energy_eV_0 = energy_eV_0
        self.sey_parameter_0 = sey_parameter_0
        self.work_function = work_function
        self.R0 = R0
        self.sey_file= sey_file
        self.energy_eV_min = energy_eV.min()
        self.energy_eV_max = energy_eV.max()
        self.delta_e = delta_e
        self.sey_diff = np.append(np.diff(sey_parameter), 0.)

    def SEY_process(self,nel_impact,E_impact_eV, costheta_impact, i_impact):

        delta = np.zeros_like(E_impact_eV)

        ref_frac = (E_impact_eV < self.work_function)
        mask_fit = (E_impact_eV > self.energy_eV_max)
        mask_regular = ~ref_frac & ~mask_fit

        delta[ref_frac] = self.R0
        delta[mask_regular] = self.interp(E_impact_eV[mask_regular])
        delta[mask_fit] = self.extrapolate_const + self.extrapolate_grad * E_impact_eV[mask_fit]

        return delta, ref_frac

    def interp(self, energy_eV):
        index_float = (energy_eV - self.energy_eV_min)/self.delta_e
        index_remainder, index_int = np.modf(index_float)
        return self.sey_parameter[index_int] + index_remainder*self.sey_diff[index_int]

    def interp_regular(self, energy_eV):
        """
        This fails if the input is not in ascending order!
        """
        print('Warning! Do not use interp_regular!')
        return np.interp(energy_eV, self.energy_eV, self.sey_parameter)

# ideas for reflected photons:


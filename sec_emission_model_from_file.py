from __future__ import division
import os
import numpy as np
import scipy.interpolate as interpolate

class SEY_from_file(object):

    def __init__(self, sey_file, R0, default_work_function=4.25):

        self.R0 = R0
        self.sey_file= sey_file

        work_function = default_work_function
        energy_eV = []
        sey_parameter = []
        with open(os.path.expanduser(sey_file)) as f:
            for ctr, line in enumerate(f):
                split = line.split()
                if 'Work function' in line:
                    work_function = float(split[-1])
                elif len(split) == 2:
                    energy_eV.append(split[0])
                    sey_parameter.append(split[1])

        self.energy_eV = np.array(energy_eV, dtype=float) + work_function
        self.sey_parameter = np.array(sey_parameter, dtype=float)
        self.work_function = work_function

        self.get_sey = interpolate.interp1d(self.energy_eV, self.sey_parameter)

#    def get_sey(self, energy_eV):
#        return np.interp(energy_eV, self.energy_eV, self.sey_parameter)

    def SEY_process(self,nel_impact,E_impact_eV, costheta_impact, i_impact):

        #E0 = self.E0
        #sq_e= np.sqrt(E_impact_eV)
        #sq_ee0 = np.sqrt(E_impact_eV+E0)
        #reflected=R0*((sq_e-sq_ee0)/(sq_e+sq_ee0))**2.

        ref_frac = (E_impact_eV < self.work_function)

        delta = np.zeros_like(E_impact_eV)
        delta[~ref_frac] = self.get_sey(E_impact_eV[~ref_frac])
        delta[ref_frac] = self.R0


        return delta, ref_frac




#
#if __name__ == '__main__':
#    obj = SEY_from_file('~/pyecloud/sey_measurements/SEY-LE_SEY/LE_SEY_Beam_Screen_V2_sample_17cm_Top.txt')


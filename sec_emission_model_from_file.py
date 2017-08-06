import numpy as np
import os

class SEY_from_file(object):

    def __init__(self, file_):

        work_function = 0.
        energy_eV = []
        sey_parameter = []
        with open(os.path.expanduser(file_)) as f:
            for ctr, line in enumerate(f):
                split = line.split()
                if 'Work function' in line:
                    work_function = float(split[-1])
                elif len(split) == 2:
                    energy_eV.append(split[0])
                    sey_parameter.append(split[1])

        self.energy_eV = np.array(energy_eV, dtype=float) + work_function
        self.sey_parameter = np.array(sey_parameter, dtype=float)

    def get_sey(self, energy_eV):
        return np.interp(energy_eV, self.energy_eV, self.sey_parameter)




if __name__ == '__main__':
    obj = SEY_from_file('~/pyecloud/sey_measurements/SEY-LE_SEY/LE_SEY_Beam_Screen_V2_sample_17cm_Top.txt')

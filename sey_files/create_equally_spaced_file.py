from __future__ import division, print_function
import numpy as np
import scipy.io as sio

def main(sey_file, range_extrapolate_right, delta_e, max_sey, output_file):
    """
    - sey file is the path to a file of the correct format, either an absolute path or in the sey_files/SEY-LE_SEY folder.
        See the example files in sey_files/SEY-LE_SEY for the format that should be used for input files.
    - range_extrapolate_right states the range in electron volts which is used to create a linear fit in order to extrapolate high energies.
    - delta_e is the resolution used for the interpolation, for example .1 eV.
    - max_sey is a factor applied to the SEY from the file, for example to emulate scrubbing effects.
        Set to None or False to disable.
    - output_file requires no explanation
    """

    # Read original file
    energy_eV_list = []
    sey_parameter_list = []
    with open(sey_file) as f:
        for ctr, line in enumerate(f):
            split = line.split()
            if line.startswith('#'):
                continue
            elif len(split) == 2:
                try:
                    energy_eV_list.append(float(split[0]))
                    sey_parameter_list.append(float(split[1]))
                except:
                    print('Error in line %i of file %s: %s' % (ctr, sey_file, line))
                    raise

    energy_eV_0 = np.array(energy_eV_list, dtype=float)
    sey_parameter_0 = np.array(sey_parameter_list, dtype=float)
    if max_sey:
        index_min_sey = np.argmin(sey_parameter_0)
        sey_parameter_0[index_min_sey:] *= max_sey/sey_parameter_0.max()

    # Build equally spaced arrays that are used by the interp function
    energy_eV = np.arange(energy_eV_0.min(), energy_eV_0.max()+delta_e*.5, delta_e)
    sey_parameter = np.interp(energy_eV, energy_eV_0, sey_parameter_0)

    # Find fit parameters for energies that extend the measured ranges
    mask_xx_fit = energy_eV_0 > (energy_eV.max() - range_extrapolate_right)
    xx_fit = energy_eV_0[mask_xx_fit]
    if len(xx_fit) < 2:
        raise ValueError('Range for SEY fit is too small! You may have to increase the range_extrapolate_right parameter!')
    yy_fit = sey_parameter_0[mask_xx_fit]
    extrapolate_grad, extrapolate_const = np.polyfit(xx_fit, yy_fit, 1)

    # save file
    save_dict = {
        'extrapolate_grad': extrapolate_grad,
        'extrapolate_const': extrapolate_const,
        'energy_eV': energy_eV,
        'sey_parameter': sey_parameter,
    }
    sio.savemat(output_file, save_dict)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('--delta', type=float, help='Energy step', default=0.1)
    parser.add_argument('--range-extrapolate-right', type=float, help='Range, starting from the highest energy, to make a linear fit for extrapolating', default=300)
    parser.add_argument('--max-sey', type=float, help='Set the maximum SEY. Values larger than the minimum of the curve are shifted accordingly')
    args = parser.parse_args()

    main(args.input, args.range_extrapolate_right, args.delta, args.max_sey, args.output)


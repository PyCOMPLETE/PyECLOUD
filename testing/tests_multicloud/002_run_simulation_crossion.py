
import sys
import os
import time
import argparse
BIN = os.path.expanduser("../../../")  # folder containing PyECLOUD, PyPIC, PyKLU
if BIN not in sys.path:
    sys.path.append(BIN)

from PyECLOUD.buildup_simulation import BuildupSimulation

sim_folder = 'LHC_ArcDriftReal_450GeV_electron_n2_sey1.75_P3.125e+06'

# check if user provided folder as command line argument
parser = argparse.ArgumentParser()
parser.add_argument('--folder', help='Simulation_folder')
parser.add_argument('--angle-dist-func',
                    help='Angular distribution of new MPs relative to surface normal. Introduced in July 2017.',
                    choices=('2D', '3D'), default='3D')

args = parser.parse_args()
if args.folder:
    sim_folder = args.folder

angle_distribution = 'cosine_%s' % args.angle_dist_func
filen_main_outp = sim_folder + '/Pyecltest_angle%s.mat' % args.angle_dist_func


time_0 = time.time()
sim = BuildupSimulation(pyecl_input_folder=sim_folder, filen_main_outp=filen_main_outp,
                        secondary_angle_distribution=angle_distribution, photoelectron_angle_distribution=angle_distribution)
sim.run()

time_needed = time.time() - time_0


print('')
print('Test simulation done in %.2f s!' % time_needed)
print('To inspect the results you can run:')
print('001_comparison_against_reference.py')
print('')

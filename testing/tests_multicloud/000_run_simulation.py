
import sys
import os
import time
import argparse

from PyECLOUD.buildup_simulation import BuildupSimulation

sim_folder = 'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_multispecies'

# check if user provided folder as command line argument
parser = argparse.ArgumentParser()
parser.add_argument('--folder', help='Simulation_folder')

args = parser.parse_args()
if args.folder:
    sim_folder = args.folder

angle_distribution = 'cosine_3D'
filen_main_outp = sim_folder + '/Pyecltest_angle3D.mat'


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

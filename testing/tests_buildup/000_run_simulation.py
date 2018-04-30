from __future__ import division, print_function
import sys
import os
import time
import argparse
BIN = os.path.expanduser("../../../") #folder containing PyECLOUD, PyPIC, PyKLU
if BIN not in sys.path:
    sys.path.append(BIN)

from PyECLOUD.buildup_simulation import BuildupSimulation

sim_folder = 'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns'
#sim_folder = 'LHC_ArcDipReal_450GeV_sey1.00_2.5e11ppb_bl_1.00ns_gas_ionization'
#sim_folder = 'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_change_s_and_E0'
#sim_folder = 'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_multigrid'
#sim_folder = 'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns'
#sim_folder = 'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_circular'
#sim_folder = 'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_skew'
#sim_folder = 'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_skew_circular'
#sim_folder = 'LHC_Sextupole_450GeV_sey1.65_2.5e11ppb_bl_1.00ns'
#sim_folder = 'LHC_Sextupole_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_skew'
#sim_folder = 'LHC_ArcDipReal_6500GeV_sey_1.70_1.1e11ppb_b1_1.00ns'
#sim_folder = 'LHC_Drift_6500GeV_sey_1.70_1.1e11ppb_b1_1.00ns'
#sim_folder = 'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_seyfromfile'
#sim_folder = 'LHC_Triplet_Quadrupole_multiple_beams'
#sim_folder = 'LHC_TDIS_non_unif_sey'
#sim_folder = 'LHC_Solenoid_sey1.10_100.00mT'
#sim_folder = 'CLIC_DRe-_Drift_0.5ns_4.0e9ppb_gas_ionization_ecloud_sey2.0'
#sim_folder = 'CLIC_DRe+_Drift_0.5ns_4.0e9ppb_gas_ionization_ecloud_sey2.0'
#sim_folder = 'CLIC_DRe-_Drift_0.5ns_4.0e9ppb_gas_ionization_ions_A18'
#sim_folder = 'CLIC_DRe+_Drift_0.5ns_4.0e9ppb_gas_ionization_ions_A18'
#sim_folder = 'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_stress_saver'

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
filen_main_outp = sim_folder+'/Pyecltest_angle%s.mat' % args.angle_dist_func


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

import sys, os
import time
os.chdir(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.expanduser("../../../") #folder containing PyECLOUD, PyPIC, PyKLU
sys.path.append(BIN)
import argparse

#ignore warnings
#import warnings
#warnings.filterwarnings("ignore")

all_sim_folders = [
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns',
    'LHC_ArcDipReal_450GeV_sey1.00_2.5e11ppb_bl_1.00ns_gas_ionization',
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_multigrid',
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_change_s_and_E0',
    'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns',
    'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_skew',
    'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_skew_circular',
    'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_circular',
    ]

parser = argparse.ArgumentParser()
parser.add_argument('index', help='Index for %s' % all_sim_folders, type=int)
args = parser.parse_args()
sim_folder = all_sim_folders[args.index]


from PyECLOUD.buildup_simulation import BuildupSimulation


time_0 = time.time()
sim = BuildupSimulation(pyecl_input_folder = sim_folder, filen_main_outp = sim_folder+'/Pyecltest.mat')
sim.run()

time_needed = time.time() - time_0


print ''
print 'Test simulation done in %.2f s!' % time_needed
print 'To inspect the results you can run:'
print '001_comparison_against_reference.py'
print ''

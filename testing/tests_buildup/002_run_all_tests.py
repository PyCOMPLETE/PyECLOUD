import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--angle-dist-func', help='Angular distribution of new MPs relative to surface normal. Introduced in July 2017.', choices=('2D', '3D'), default='3D')
parser.add_argument('--all', help='Run 2D and 3D consecutively', choices=('2D', '3D'), default='3D')
args = parser.parse_args()


all_sim_folders = [
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns',
    'LHC_ArcDipReal_450GeV_sey1.00_2.5e11ppb_bl_1.00ns_gas_ionization',
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_multigrid',
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_change_s_and_E0',
    'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns',
    'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_skew',
    'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_circular',
    'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_skew_circular',
    'LHC_Sextupole_450GeV_sey1.65_2.5e11ppb_bl_1.00ns',
    'LHC_Sextupole_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_skew',
    ]


if args.all:
    dists = ['2D', '3D']
else:
    dists = [args.angle_dist_func]

for ctr, sim_folder in enumerate(all_sim_folders):
    for dist in dists:
        for cmd in [
            'python2 ./000_run_simulation.py --folder %s --angle-dist-func %s' % (sim_folder, dist),
            'python2 ./001_comparison_against_reference.py --folder %s --angle-dist-func %s' % (sim_folder, dist)
        ]:
            status = os.system(cmd)
            if status != 0:
                print('%s finished with status %i' % (cmd, status))
                sys.exit()


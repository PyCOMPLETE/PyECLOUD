import os
import argparse

all_dists = ('2D', '3D')

parser = argparse.ArgumentParser()
parser.add_argument('--angle-dist-func', help='Angular distribution of new MPs relative to surface normal. Introduced in July 2017.', choices=all_dists, default='3D')
parser.add_argument('--all', help='Run 2D and 3D consecutively', action='store_true')
args = parser.parse_args()

all_sim_folders = [
    'LHC_ArcDipReal_450GeV_sey1.00_2.5e11ppb_bl_1.00ns_gas_ionization',
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns',
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_change_s_and_E0',
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_seyfromfile',
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_stress_saver',
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_multigrid',
    'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns',
    'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_circular',
    'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_skew',
    'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_skew_circular',
    'LHC_Sextupole_450GeV_sey1.65_2.5e11ppb_bl_1.00ns',
    'LHC_Sextupole_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_skew',
    'LHC_Drift_6500GeV_sey_1.70_1.1e11ppb_b1_1.00ns',
    'LHC_ArcDipReal_6500GeV_sey_1.70_1.1e11ppb_b1_1.00ns',
    'LHC_Triplet_Quadrupole_multiple_beams',
    'CLIC_DRe-_Drift_0.5ns_4.0e9ppb_gas_ionization_ecloud_sey2.0',
    'CLIC_DRe+_Drift_0.5ns_4.0e9ppb_gas_ionization_ecloud_sey2.0',
    'CLIC_DRe-_Drift_0.5ns_4.0e9ppb_gas_ionization_ions_A18',
    'CLIC_DRe+_Drift_0.5ns_4.0e9ppb_gas_ionization_ions_A18',
    'LHC_TDIS_non_unif_sey',
    'LHC_Solenoid_sey1.10_100.00mT',
    'FCC_Dipole_25ns_50.00TeV_sey1.9_segment_photoemission',
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_nonuniftime',
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_checkpoint',
    'Rectangular_Dip_450GeV_sey1.60_1.1e11ppb_furman_pivi',
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_em_tracking'
    ]

if args.all:
    dists = all_dists
else:
    dists = [args.angle_dist_func]

for ctr, sim_folder in enumerate(all_sim_folders):
    test_script = './000_run_simulation.py'
    if sim_folder.endswith('checkpoint'):
        test_script = './000b_run_test_checkpoint.py'
    for dist in dists:
        for cmd in [
            'python %s --folder %s --angle-dist-func %s' % (test_script, sim_folder, dist),
            'python ./001_comparison_against_reference.py --folder %s --angle-dist-func %s' % (sim_folder, dist)
        ]:
            print(cmd)
            status = os.system(cmd)
            if status != 0:
                raise SystemError('%s finished with status %i' % (cmd, status))

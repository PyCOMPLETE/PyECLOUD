import shutil
import os

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
    'LHC_Octupole_6500GeV_sey1.65_2.5e11ppb_b1_1.00ns',
    'LHC_Octupole_6500GeV_sey1.65_2.5e11ppb_b1_1.00ns_skew',
    ]


for folder in all_sim_folders:
    ref_file = folder + '/Pyecltest_ref.mat'
    mat_file = folder + '/Pyecltest.mat'
    os.remove(ref_file)
    shutil.copy(mat_file, ref_file)
    print('%s replaced by %s' % (ref_file, mat_file))


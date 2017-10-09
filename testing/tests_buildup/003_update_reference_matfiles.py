import sys
import shutil
import os

if sys.version_info.major != 2:
    raw_input = input

if raw_input('Continue replace all reference files? y/n ') not in ('y','yes'):
    sys.exit()

all_sim_folders = [
    'LHC_ArcDipReal_450GeV_sey1.00_2.5e11ppb_bl_1.00ns_gas_ionization',
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns',
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_change_s_and_E0',
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
]

for folder in all_sim_folders:
    for dim in 2, 3:
        ref_file = folder + '/Pyecltest_angle%iD_ref.mat' % dim
        mat_file = folder + '/Pyecltest_angle%iD.mat' % dim
        if os.path.isfile(ref_file):
            os.remove(ref_file)
        shutil.copy(mat_file, ref_file)
        print('%s replaced by %s' % (ref_file, mat_file))

# OLD STUFF:
#for folder in all_sim_folders:
#    ref_file = folder + '/Pyecltest_ref.mat'
#    new_ref_file = folder + '/Pyecltest_3D_ref.mat'
#    ref_file_2D = './old_tests_cosine_2D/' + ref_file
#    new_ref_file_2D = folder + 'Pyecltest_2D_ref.mat'
#    shutil.move(ref_file, new_ref_file)
#    shutil.copy(ref_file_2D, new_ref_file_2D)
#    print('%s replaced by %s' % (ref_file, new_ref_file))
#    print('%s replaced by %s' % (ref_file_2D, new_ref_file_2D))

#for folder in all_sim_folders:
#    new_ref_file = folder + 'Pyecltest_3D_ref.mat'
#    ref_file = folder + '/Pyecltest_3D_ref.mat'
#    new_ref_file_2D = folder + 'Pyecltest_2D_ref.mat'
#    ref_file_2D = folder + '/Pyecltest_2D_ref.mat'
#    shutil.move(new_ref_file, ref_file)
#    shutil.move(new_ref_file_2D, ref_file_2D)

#for folder in all_sim_folders:
#    old_mat = folder +'/Pyecltest_3D_ref.mat'
#    new_mat = folder +'/Pyecltest_angle3D_ref.mat'
#    shutil.move(old_mat, new_mat)

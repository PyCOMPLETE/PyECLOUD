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
    for dim in 2, 3:
        ref_file = folder + '/Pyecltest_angle_%iD_ref.mat' % dim
        mat_file = folder + '/Pyecltest_angle_%iD.mat' % dim
        os.remove(ref_file)
        shutil.copy(mat_file, ref_file)
        print('%s replaced by %s' % (ref_file, mat_file))


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

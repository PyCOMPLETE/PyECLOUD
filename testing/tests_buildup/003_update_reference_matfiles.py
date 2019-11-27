import sys
import shutil
import os

if sys.version_info.major != 2:
    raw_input = input

if input('Continue replace all reference files? y/n ') not in ('y', 'yes'):
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
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_stress_saver'
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
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_seyfromfile',
    'LHC_Solenoid_sey1.10_100.00mT',
    'FCC_Dipole_25ns_50.00TeV_sey1.9_segment_photoemission',
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_nonuniftime'
]

for folder in all_sim_folders:
    #for dim in 2, 3:
    for dim in [3]:
        ref_file = folder + '/Pyecltest_angle%iD_ref.mat' % dim
        mat_file = folder + '/Pyecltest_angle%iD.mat' % dim
        if os.path.isfile(ref_file):
            os.remove(ref_file)
        shutil.copy(mat_file, ref_file)
        print(('%s replaced by %s' % (ref_file, mat_file)))

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

#dist = 'cosine_3D'
#for sim in all_sim_folders:
#    filename = 'secondary_emission_parameters.input'
#    with open(os.path.join(sim, filename), 'a') as f:
#        f.write("\n\nsecondary_angle_distribution = '%s'\n" % dist)
#    print('Modified %s in %s' % (filename, sim))
#    filename = 'machine_parameters.input'
#    with open(os.path.join(sim, filename), 'a') as f:
#        f.write("\n\nphotoelectron_angle_distribution = '%s'\n" % dist)
#    print('Modified %s in %s' % (filename, sim))
#

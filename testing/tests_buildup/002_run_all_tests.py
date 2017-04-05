import os

all_sim_folders = [
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns',
    'LHC_ArcDipReal_450GeV_sey1.00_2.5e11ppb_bl_1.00ns_gas_ionization',
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_multigrid',
    'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_change_s_and_E0',
    'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns',
    ]



for ctr, sim_folder in enumerate(all_sim_folders):
    os.system('python ./000_run_simulation.py %i' % ctr)
    os.system('rm -r %s/comparison_plots' % sim_folder)
    os.system('python ./001_comparison_against_reference.py %i' % ctr)



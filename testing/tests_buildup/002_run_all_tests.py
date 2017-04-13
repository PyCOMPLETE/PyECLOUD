import os
from tests import all_sim_folders



for ctr, sim_folder in enumerate(all_sim_folders):
    os.system('python ./000_run_simulation.py %i' % ctr)
    os.system('rm -r %s/comparison_plots' % sim_folder)
    os.system('python ./001_comparison_against_reference.py %i' % ctr)



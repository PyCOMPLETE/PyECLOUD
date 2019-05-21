import numpy as np

from PyECLOUD.buildup_simulation import BuildupSimulation

sim_input_folder = './'
t_stop = None

sim = BuildupSimulation(pyecl_input_folder=sim_input_folder, 
        filen_main_outp='./Pyecltest.mat',
        extract_sey=False, 
        )

ec = sim.cloud_list[0]

sim.run(t_end_sim = t_stop)
    

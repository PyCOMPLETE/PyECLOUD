import numpy as np

from PyECLOUD.buildup_simulation import BuildupSimulation

sim_input_folder = './'
t_stop = None

save_once_custom_observables = {
        'Emax': lambda sim: sim.cloud_list[0].impact_man.sey_mod.Emax,
        'del_max': lambda sim: sim.cloud_list[0].impact_man.sey_mod.del_max,
        'R0': lambda sim: sim.cloud_list[0].impact_man.sey_mod.R0,
        'E0': lambda sim: sim.cloud_list[0].impact_man.sey_mod.E0,
        's': lambda sim: sim.cloud_list[0].impact_man.sey_mod.s,
        'sigmafit': lambda sim: sim.cloud_list[0].impact_man.sey_mod.sigmafit,
        'mufit': lambda sim: sim.cloud_list[0].impact_man.sey_mod.mufit
         }

sim = BuildupSimulation(pyecl_input_folder=sim_input_folder, 
        filen_main_outp='./Pyecltest.mat',
        extract_sey=False,
        save_once_custom_observables = save_once_custom_observables
        )

ec = sim.cloud_list[0]

sim.run(t_end_sim = t_stop)
    

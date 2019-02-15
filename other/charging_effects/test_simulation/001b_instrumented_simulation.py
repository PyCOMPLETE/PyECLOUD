from PyECLOUD.buildup_simulation import BuildupSimulation

import numpy as np

step_by_step_custom_observables = {
        'Qpatch_ave': lambda sim: np.mean(sim.cloud_list[0].impact_man.sey_mod.Q_segments)
        }

sim = BuildupSimulation(
        step_by_step_custom_observables=step_by_step_custom_observables)

sim.run(t_end_sim = None) 


from PyECLOUD.buildup_simulation import BuildupSimulation

import numpy as np

step_by_step_custom_observables = {
        'Qpatch_ave': lambda sim: np.mean(sim.cloud_list[0].impact_man.sey_mod.Q_segments[
            sim.cloud_list[0].impact_man.sey_mod.flag_charging])
        }
pass_by_pass_custom_observables = {
        'Q_segments' : lambda sim: sim.cloud_list[0].impact_man.sey_mod.Q_segments.copy()
        }

sim = BuildupSimulation(
        step_by_step_custom_observables=step_by_step_custom_observables,
        pass_by_pass_custom_observables=pass_by_pass_custom_observables)

sim.run(t_end_sim = None) 


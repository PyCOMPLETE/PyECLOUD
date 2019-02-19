from PyECLOUD.buildup_simulation import BuildupSimulation

import numpy as np

def sey_at_emax_patch(sim):
    ec = ec = sim.cloud_list[0]
    flag_patch = ec.impact_man.sey_mod.flag_charging
    i_patch = np.where(flag_patch)[0]
    Emax_patch = ec.impact_man.sey_mod.Emax_segments[flag_patch]
    
    nel_probe = 0.0001
    nel_out, _, _ = ec.impact_man.sey_mod.SEY_process(nel_impact=0*Emax_patch+nel_probe, 
            E_impact_eV=Emax_patch, 
            costheta_impact=0*Emax_patch+1.,
            i_impact = i_patch)
    del_emax = np.mean(nel_out)/nel_probe

    return del_emax


step_by_step_custom_observables = {
        'Qpatch_ave': lambda sim: np.mean(sim.cloud_list[0].impact_man.sey_mod.Q_segments[
            sim.cloud_list[0].impact_man.sey_mod.flag_charging]),
        'sey_at_emax_patch': sey_at_emax_patch
        }
pass_by_pass_custom_observables = {
        'Q_segments' : lambda sim: sim.cloud_list[0].impact_man.sey_mod.Q_segments.copy()
        }

sim = BuildupSimulation(
        step_by_step_custom_observables=step_by_step_custom_observables,
        pass_by_pass_custom_observables=pass_by_pass_custom_observables)

sim.run(t_end_sim = None) 


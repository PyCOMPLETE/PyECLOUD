import numpy as np

from PyECLOUD.buildup_simulation import BuildupSimulation

sim_input_folder = '../../testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns'

t_stop_list = [5e-9, 25e-9, 123e-9, 25e-9, 250e-9, 500e-9, None]


def observe_Nel(sim):
    N_mp = sim.cloud_list[0].MP_e.N_mp
    Nel = np.sum(sim.cloud_list[0].MP_e.nel_mp[:N_mp])
    return Nel

step_by_step_custom_observables = {
        'Nelectrons': observe_Nel,
        'y_first_ele': lambda sim: sim.cloud_list[0].MP_e.y_mp[0],
        'relative_charge_first_ele': lambda sim: sim.cloud_list[0].MP_e.nel_mp[0]/sim.cloud_list[0].MP_e.nel_mp_ref,
        'ref_size': lambda sim: sim.cloud_list[0].MP_e.nel_mp_ref
        }

pass_by_pass_custom_observables = {
        'sum_rho': lambda sim: np.sum(sim.cloud_list[0].rho, axis=1)
        }

save_once_custom_observables = {
        'Vx': lambda sim: sim.cloud_list[0].impact_man.chamb.Vx
        }

sim = BuildupSimulation(pyecl_input_folder=sim_input_folder, 
        filen_main_outp='./Pyecltest.mat',
        extract_sey=False, 
        step_by_step_custom_observables = step_by_step_custom_observables,
        pass_by_pass_custom_observables = pass_by_pass_custom_observables, 
        save_once_custom_observables = save_once_custom_observables)

ec = sim.cloud_list[0]

for ii, t_stop in enumerate(t_stop_list):
    print('\n\n==============================')
    print(('Simulation run %d - t_stop = %s s'%(ii, repr(t_stop)))) 
    print((' starting at tt=%s s'%repr(sim.beamtim.tt_curr))) 
    
    sim.run(t_end_sim = t_stop)
    
    print((' after run a tt=%s'%repr(sim.beamtim.tt_curr)))

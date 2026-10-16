import numpy as np

from PyECLOUD.buildup_simulation import BuildupSimulation

sim_input_folder = '../../testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns'
sim_input_folder = 'input_files'
t_stop_list = [1800e-9]

i_observe = range(0, 60000, 100)

def observe_electrons(sim):

    ec = sim.cloud_list[0]

    if not hasattr(ec, 'x_follow'):
        ec.x_follow = []
        ec.y_follow = []
        ec.z_follow = []
        ec.nel_follow = []

    mask_observe = np.zeros_like(ec.MP_e.x_mp, dtype=np.bool)
    mask_observe[i_observe] = True

    x_temp = ec.MP_e.x_mp[mask_observe]
    y_temp = ec.MP_e.y_mp[mask_observe]
    z_temp = ec.MP_e.z_mp[mask_observe]
    nel_temp = ec.MP_e.nel_mp[mask_observe]

    temp_mask = mask_observe.copy()
    temp_mask[ec.MP_e.N_mp:] = False
    N_keep = np.sum(temp_mask)

    x_temp[N_keep:] = np.nan
    y_temp[N_keep:] = np.nan
    z_temp[N_keep:] = np.nan
    nel_temp [N_keep:] = np.nan

    ec.x_follow.append(x_temp)
    ec.y_follow.append(y_temp)
    ec.z_follow.append(z_temp)
    ec.nel_follow.append(nel_temp)

    return np.sum(nel_temp)

step_by_step_custom_observables = {
        'N_observed': observe_electrons,
        }

pass_by_pass_custom_observables = {
        }

save_once_custom_observables = {
        }

sim = BuildupSimulation(pyecl_input_folder=sim_input_folder,
        logfile_path='logfile.txt',
        progress_path='progress.txt',
        filen_main_outp='./Pyecltest.mat',
        secondary_angle_distribution = 'cosine_3D',
        fact_clean = 0.,  # cleanings would move electrons around
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


x_foll = np.array(ec.x_follow)
y_foll = np.array(ec.y_follow)
z_foll = np.array(ec.z_follow)
nel_foll = np.array(ec.nel_follow)

import scipy.io as sio
sio.savemat('followed_electrons.mat',{
    'x': x_foll,
    'y': y_foll,
    'z': z_foll,
    'nel': nel_foll,
    }, oned_as='row')

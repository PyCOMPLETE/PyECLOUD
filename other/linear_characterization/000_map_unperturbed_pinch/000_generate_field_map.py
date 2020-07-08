import os

import numpy as np
import scipy.io as sio

from PyPARIS_sim_class import Simulation as sim_mod
import PyPARIS.util as pu

sim_param_file = '../reference_simulation/Simulation_parameters.py'
sim_param_amend_files = ['../Simulation_parameters_amend.py']

# Instantiate simulation
sim_content = sim_mod.Simulation(param_file=sim_param_file)

# Here sim_content.pp can be edited (directly and through files)
for ff in sim_param_amend_files:
    sim_content.pp.update(param_file=ff)

# Add ring of CPU information (mimicking the master core)
pu.get_sim_instance(sim_content,
        N_cores_pretend=sim_content.pp.n_segments,
        id_pretend=sim_content.pp.n_segments-1,
        init_sim_objects_auto=False)
assert(sim_content.ring_of_CPUs.I_am_the_master)

# Initialize machine elements
sim_content.init_all(generate_parent_eclouds=True,
            install_clouds=False)

# Initialize slicer
if os.path.exists('simulation_status.sta'):
    os.remove('simulation_status.sta')
sim_content.init_master(generate_bunch=False,
        prepare_monitors=False)

# Generate bunch for map
pp = sim_content.pp
bunch_for_map = sim_content.machine.generate_6D_Gaussian_bunch_matched(
    n_macroparticles=pp.n_macroparticles_for_footprint_map,
    intensity=pp.intensity,
    epsn_x=pp.epsn_x,
    epsn_y=pp.epsn_y,
    sigma_z=pp.sigma_z,
)

# Slice the bunch
slicer = sim_content.slicer
slices_list_for_map = bunch_for_map.extract_slices(slicer)

# Record field maps
Ex_L_map = 0.
Ey_L_map = 0.
for iee, ee in enumerate(sim_content.parent_eclouds):
    ee.save_ele_distributions_last_track = True
    ee.track_once_and_replace_with_recorded_field_map(
        slices_list_for_map)
    Ex_L_map = Ex_L_map + ee.efieldmap.Ex * ee.efieldmap.L_interaction * sim_content.n_segments
    Ey_L_map = Ey_L_map + ee.efieldmap.Ey * ee.efieldmap.L_interaction * sim_content.n_segments
    sio.savemat(f'rho_map_ec{iee}.mat', {
        'xg': ee.efieldmap.pic.xg,
        'yg': ee.efieldmap.pic.yg,
        'zg': [ss.slice_info['z_bin_center'] for ss in slices_list_for_map],
        'rho': ee.rho_ele_last_track,
        'sigma_x': bunch_for_map.sigma_x(),
        'sigma_y': bunch_for_map.sigma_y()
        })

sio.savemat('field_map.mat', {
    'xg': ee.efieldmap.pic.xg,
    'yg': ee.efieldmap.pic.yg,
    'zg': [ss.slice_info['z_bin_center'] for ss in slices_list_for_map],
    'Ex_L_map': Ex_L_map,
    'Ey_L_map': Ey_L_map,
    'sigma_x': bunch_for_map.sigma_x(),
    'sigma_y': bunch_for_map.sigma_y()
    })

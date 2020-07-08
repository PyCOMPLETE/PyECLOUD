from scipy.constants import c

# check_for_resubmit = False

####################
# Machine Settings #
####################

machine_configuration = 'HLLHC-injection'

# # Use this part for optics from file
# # n_segments needs to be None if optics_pickle_file is specified
# optics_pickle_file = 'lhc2018_25cm_only_triplets_IR15_b1_optics.pkl'
# n_segments = None
# beta_x =  None
# beta_y =  None
# Q_x = None
# Q_y = None

# # Use this part for smooth machine
optics_pickle_file = None
n_segments = 8
beta_x =  92.7
beta_y =  93.2
Q_x = 62.27
Q_y = 60.295

Qp_x = 0.
Qp_y = 0.

octupole_knob = 0.

n_non_parallelizable = 2 #rf and aperture

V_RF = 6e6


# Transverse Damper Settings
enable_transverse_damper = False
dampingrate_x = 100.
dampingrate_y = 100.
if enable_transverse_damper: n_non_parallelizable += 1


###################
# Beam Parameters #
###################

bunch_from_file = None

intensity = 1.2e+11

epsn_x = 2.5e-6
epsn_y = 2.5e-6

sigma_z = 0.08656507 #1.2e-9/4*c

x_kick_in_sigmas = 0.1
y_kick_in_sigmas = 0.1

# Numerical Parameters
n_slices = 500
z_cut = 2.5e-9/2*c # For slicing
macroparticles_per_slice = 5000
n_macroparticles = macroparticles_per_slice*n_slices


#################
# Stop Criteria #  
#################

# 1. Turns
N_turns = 128 # Per job
N_turns_target = 20000
# 2. Losses
sim_stop_frac = 0.9
# 3. Emittance Growth
flag_check_emittance_growth = True
epsn_x_max_growth_fraction = 0.5
epsn_y_max_growth_fraction = epsn_x_max_growth_fraction


######################
# Footprint Settings #
######################

footprint_mode = False
n_macroparticles_for_footprint_map = 500000
n_macroparticles_for_footprint_track = 5000


####################
# E-Cloud Settings #
####################

# General E-Cloud Settings
chamb_type = 'polyg'
x_aper = 2.300000e-02
y_aper = 1.800000e-02
filename_chm = 'LHC_chm_ver.mat'
Dt_ref = 5e-12
pyecl_input_folder = './pyecloud_config'
sey = 1.40

# Transverse Multigrid Parameters
PyPICmode = 'ShortleyWeller_WithTelescopicGrids'
N_min_Dh_main = 10.
Dh_sc_ext = .8e-3
f_telescope = 0.3
N_nodes_discard = 5.
target_size_internal_grid_sigma = 10.
target_Dh_internal_grid_sigma = 0.2
custom_target_grid_arcs = None

# # Uncomment for custom grid
# custom_target_grid_arcs = {
#     'x_min_target': -3e-3,
#     'x_max_target': 3e-3,
#     'y_min_target': -3.1e-3,
#     'y_max_target': 3.1e-3,
#     'Dh_target': 7e-5}

force_interp_at_substeps_interacting_slices = True

# Enable Kicks Different Planes
enable_kick_x = True
enable_kick_y = True

# Dedicated Dipole E-Cloud Settings
enable_arc_dip = True
fraction_device_dip = 6.000000e-01
init_unif_edens_flag_dip = 1
init_unif_edens_dip = 12.000000e+11
N_MP_ele_init_dip = 1000000
N_mp_max_dip = N_MP_ele_init_dip*4
B_multip_dip = [0.54] #[8.33] #T

# Dedicated Quadrupole E-Cloud Settings
enable_arc_quad = False
fraction_device_quad = 7.000000e-02
N_mp_max_quad = 2000000 
B_multip_quad = [0., 12.1] #T
folder_path = '../../LHC_ecloud_distrib_quads/'
eMPs = 500000
filename_state = 'combined_distribution_sey_%.2f_sigmat_%.3fns_450Gev_N_mp_%d_symm'%(sey, sigma_z/c*1e9,eMPs)
filename_init_MP_state_quad = folder_path + filename_state

# Dedicated Kick Element Settings
enable_eclouds_at_kick_elements = False
path_buildup_simulations_kick_elements = '/home/kparasch/workspace/Triplets/ec_headtail_triplets/simulations_PyECLOUD/!!!NAME!!!_sey1.35'
name_MP_state_file_kick_elements = 'MP_state_9.mat'
orbit_factor = 6.250000e-01


from scipy.constants import c
from scipy.constants import m_p
from scipy.constants import e as qe

####################
# Machine Settings #
####################

machine_class = 'Synchrotron'

optics_mode = 'smooth'
charge = qe
mass = m_p
p0 = 450e9 * qe / c
circumference = 26658.8832
n_segments = 16
name = None
s = None
alpha_x = 0.
beta_x = 92.7
D_x = 0.
alpha_y = 0
beta_y = 93.2
D_y = 0.
accQ_x = 62.27
accQ_y = 60.295
Qp_x = 0.
Qp_y = 0.
app_x = 0.
app_y = 0.
app_xy = 0.
longitudinal_mode = 'non-linear'
Q_s = None
alpha_mom_compaction = 3.225e-04
h_RF = 35640
V_RF = 8.0e6
dphi_RF = 0.
p_increment = 0.
RF_at = 'end_of_transverse'
wrap_z = False
other_detuners = []

# Transverse Damper Settings
enable_transverse_damper = True
dampingrate_x = 50.
dampingrate_y = 100.

###################
# Beam Parameters #
###################

bunch_from_file = None

intensity = 1.2e+11

epsn_x = 2.5e-6
epsn_y = 2.5e-6

sigma_z = 1.2e-9/4*c

x_kick_in_sigmas = 0.1
y_kick_in_sigmas = 0.1

# Numerical Parameters
n_slices = 500
z_cut = 2.5e-9/2*c # For slicing
macroparticles_per_slice = 100
n_macroparticles = macroparticles_per_slice*n_slices


#################
# Stop Criteria #  
#################

# 1. Turns
N_turns = 50 # Per job
N_turns_target = 150 
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
sey = 1.30

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
enable_arc_dip = False
fraction_device_dip = 0.65
init_unif_edens_flag_dip = 1
init_unif_edens_dip = 1.000000e+12
N_MP_ele_init_dip = 500000
N_mp_max_dip = N_MP_ele_init_dip*4
B_multip_dip = [8.33] #T

# Dedicated Quadrupole E-Cloud Settings
enable_arc_quad = False 
fraction_device_quad = 26.000000e-02 #7.000000e-02
N_mp_max_quad = 2000000 
B_multip_quad = [0., 188.2] #T
folder_path = '../../LHC_ecloud_distrib_quads/'
filename_state = 'combined_distribution_sey%.2f_%.1fe11ppb_7tev.mat'%(sey,intensity/1e11)
filename_init_MP_state_quad = folder_path + filename_state

# Dedicated Kick Element Settings
enable_eclouds_at_kick_elements = False 
path_buildup_simulations_kick_elements = '/home/kparasch/workspace/Triplets/ec_headtail_triplets/simulations_PyECLOUD/!!!NAME!!!_sey1.35'
name_MP_state_file_kick_elements = 'MP_state_9.mat'
orbit_factor = 6.250000e-01


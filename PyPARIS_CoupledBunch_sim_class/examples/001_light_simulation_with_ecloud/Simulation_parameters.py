
resubmit_command = None

####################
# Machine Settings #
####################

machine_class_path = 'PyPARIS_CoupledBunch_sim_class.LHC_custom.LHC'
machine_configuration = 'HLLHC-injection'
n_segments = 2

Qp_x = 0.
Qp_y = 0.
octupole_knob = 0.0

# Transverse Damper Settings
enable_transverse_damper = False
dampingrate_x = 20.
dampingrate_y = 300000.


###################
# Beam Parameters #
###################

bunch_intensity = 1e11
epsn_x = 2.5e-6
epsn_y = 2.5e-6
sigma_z = 10e-2

# Beam properties
#Filling pattern: here head is left and tail is right
b_spac_s = 25e-9/5
filling_pattern = 1 * (2*([1.]+4*[0.]) )

load_beam_from_folder = None #'bunch_states_turn0'

macroparticlenumber = 100000

non_linear_long_matching = False

x_kick_in_sigmas = 0.25
y_kick_in_sigmas = 0.25


####################
# Slicing settings #
####################

n_slices_per_bunch = 100
z_cut_slicing = 3*sigma_z
min_inten_slice4EC = 1e7

#######################
# Simulation settings #
#######################

N_turns_target = 20000
N_turns = 3 
N_parellel_rings = 3
N_pieces_per_transfer = n_slices_per_bunch + 1
verbose = False
mpi_verbose = False
enable_barriers = True

###########
# Savings #
###########

save_beam_at_turns = []
stats_to_store = [
             'mean_x', 'mean_xp', 'mean_y', 'mean_yp', 'mean_z', 'mean_dp',
             'sigma_x', 'sigma_y', 'sigma_z','sigma_dp', 'epsn_x', 'epsn_y',
             'epsn_z', 'macroparticlenumber',
             'i_bunch', 'i_turn']
slice_stats_to_store='mean_x mean_y mean_z n_macroparticles_per_slice'.split()


####################
# e-cloud settings #
####################

enable_ecloud = True

enable_kick_x = True
enable_kick_y = False

L_ecloud_tot = 20e3

target_size_internal_grid_sigma = 10.
Dt_ref = 25e-12
pyecl_input_folder = './pyecloud_config'
chamb_type = 'polyg' 
filename_chm = 'LHC_chm_ver.mat'
save_only = ['lam_t_array', 'nel_hist', 'Nel_timep', 't', 't_hist', 'xg_hist']


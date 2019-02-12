from scipy.constants import m_p, m_e, e as qe
import numpy as np

# To be safe, only use immutable python types as default values in the parameters_dict.
# This most importantly excludes python lists. Tuples may be used in their place.

parameters_dict = {
    'superparameters': {'pi'},  # are allowed anywhere and can be repeated (for backwards compatibility)
    'simulation_parameters': {
        'mandatory': {

            # Other input file names
            'machine_param_file',
            'secondary_emission_parameters_file',
            'beam_parameters_file',

            # Log and progress files
            'logfile_path',
            'progress_path',

            # Time sampling
            'Dt',
            't_end',

            # Neglibile beam linear density
            'lam_th',

            # MP management settings
            'N_mp_max',
            'N_mp_regen',
            'N_mp_regen_low',
            't_ON_regen_low',
            'N_mp_after_regen',
            'fact_split',
            'fact_clean',
            'nel_mp_ref_0',
            'Nx_regen', 'Ny_regen', 'Nvx_regen', 'Nvy_regen', 'Nvz_regen',
            'regen_hist_cut',

            # Space charge parameters
            'Dt_sc',
            'Dh_sc',
            't_sc_ON',

            # Saving settings
            'Dx_hist',
            'r_center',
            'Dt_En_hist',
            'Nbin_En_hist',
            'En_hist_max',
        },
        'optional': {
            # Secondary beams
            'secondary_beams_file_list': (),

            # Additional clouds
            'additional_clouds_file_list': (),

            # Name, mass and charge for default cloud
            'cloud_name': None,
            'cloud_mass': m_e,
            'cloud_charge': -qe,

            'N_mp_soft_regen': None,
            'N_mp_after_soft_regen': None,
            'stopfile': 'stop',

            # Saving settings
            'filen_main_outp': 'Pyecltest.mat',
            'flag_movie': 0,
            'flag_sc_movie': 0,
            'flag_cos_angle_hist': True,
            'cos_angle_width' : 0.05,
            'save_mp_state_time_file': -1,
            'flag_detailed_MP_info': 0,
            'flag_hist_impact_seg': 0,
            'flag_verbose_file': False,
            'flag_verbose_stdout': False,
            'dec_fac_secbeam_prof': 1,
            'el_density_probes': (),
            'save_simulation_state_time_file': -1,
            'checkpoint_DT': None,
            'checkpoint_folder': None,
            'copy_main_outp_DT': None,
            'copy_main_outp_folder': None,
            'x_min_hist_det': None,
            'y_min_hist_det': None,
            'x_max_hist_det': None,
            'y_max_hist_det': None,
            'Dx_hist_det': None,

            'sparse_solver': 'scipy_slu',
            'PyPICmode'    : 'FiniteDifferences_ShortleyWeller',

            # Multigrid parameters
            'f_telescope': None,
            'target_grid': None,
            'N_nodes_discard': None,
            'N_min_Dh_main': None,

            'dec_fact_out': 1,

            # Where to put this?
            't_ion': -1,

            'extract_sey': True,
        },
    },
    'machine_parameters': {
        'mandatory': set(),
        'optional': {

            # Chamber profile
            'chamb_type': 'ellip',
            'x_aper': None,
            'y_aper': None,
            'filename_chm': None,
            'filename_chm_photoem': None,

            # Tracking and magnetic field
            'track_method': 'StrongBdip',
            'B': 0.,  # Tesla (if B=-1 computed from energy and bending radius)
            'bm_totlen': -1,  # m
            'B_map_file': None,
            'Bz_map_file': None, # documented?
            'fact_Bmap': 1.,
            'B0x': 0.,
            'B0y': 0.,
            'B0z': 0.,
            'B_zero_thrhld': None,
            'N_sub_steps': 1,
            'B_multip': [],
            'B_skew': None,

            # Optics
            'betafx': None,
            'betafy': None,
            'Dx': 0.,
            'Dy': 0.,

            # Residual gas ionization
            'gas_ion_flag'      : 0,
            'P_nTorr'           : -1,
            'sigma_ion_MBarn'   : -1,
            'Temp_K'            : -1,
            'unif_frac'         : -1,
            'E_init_ion'        : -1,

            # Photoemission
            'photoem_flag'                      : 0,
            'inv_CDF_refl_photoem_file'         : -1,
            'inv_CDF_all_photoem_file'          : -1,
            'k_pe_st'                           : -1,
            'refl_frac'                         : -1,
            'alimit'                            : -1,
            'e_pe_sigma'                        : -1,
            'e_pe_max'                          : -1,
            'x0_refl'                           : -1,
            'y0_refl'                           : -1,
            'out_radius'                        : -1,
            'phem_resc_fac'                     : 0.9999,
            'photoelectron_angle_distribution'  : 'undefined',
            'energy_distribution'               : 'gaussian',
            'flag_continuous_emission'          : False,


            # Uniform initial distribution
            'init_unif_flag'        : 0,
            'Nel_init_unif'         : None,
            'E_init_unif'           : 0,
            'x_max_init_unif'       : None,
            'x_min_init_unif'       : None,
            'y_max_init_unif'       : None,
            'y_min_init_unif'       : None,
            'filename_init_MP_state': None, # undocumented?

            # Uniform initial density
            'init_unif_edens_flag'  : 0,
            'init_unif_edens'       : None,
            'E_init_unif_edens'     : 0.,
            'x_max_init_unif_edens' : None,
            'x_min_init_unif_edens' : None,
            'y_max_init_unif_edens' : None,
            'y_min_init_unif_edens' : None,

            'flag_assume_convex': True,
        },
    },
    'beam_beam': {
        'mandatory': {

            # Basic definitions
            'energy_eV',

            #Transverse electric field
            'beam_field_file',

            # Beam longitudinal profile
            'b_spac',
            'fact_beam',
            'flag_bunched_beam',

            # this is mandatory!
            't_offs',
            'filling_pattern_file',
        },
        'optional': {
            # Basic definitions
            'q_part': qe,
            'm0_part': m_p,
            'Dp_p': 0.,
            'nemittx': None,
            'nemitty': None,
            'x_beam_pos': 0.,
            'y_beam_pos': 0.,
            'sigmax': -1,
            'sigmay': -1,

            #Transverse electric field
            'save_beam_field_file_as': None,

            # if beam_field_file is given
            'Dh_beam_field': None,
            'Nx': None,
            'Ny': None,
            'nimag': None,

            # if compute_FDSW_multigrid
            'Dh_beam_field': None,
            'f_telescope_beam': None,
            'target_grid_beam': None,
            'N_nodes_discard_beam': None,
            'N_min_Dh_main_beam': None,

            # if flag_bunched_beam == 1
            'sigmaz' : -1,

            # if flag_bunched_beam == 0
            'beam_long_prof_file': None,

            #????
            'Dx': None,
            'Dy': None,
            'betafx': None,
            'betafy': None,

            # this is optional!
            'coast_dens': 0.
        },
    },
    'secondary_emission_parameters': {
        'mandatory': {

            # Secondray Electron Energy Spectrum
            'E_th',
            'sigmafit',
            'mufit',

            # Other parameters
            'scrub_en_th',
        },
        'optional': {

            # Secondray Electron Yield
            'Emax': None,
            'del_max': None,
            'R0': None,
            'E0': None,
            's_param': None,

            # Choice of model
            'switch_model': 0,

            # Other parameters
            'secondary_angle_distribution': 'undefined',
            'switch_no_increase_energy': 0,
            'thresh_low_energy': -1,

            # SEY from file
            'sey_file': None,
            'flag_costheta_Emax_shift': True,
            'flag_costheta_delta_scale': True,

            # Furman-Pivi Model
            'p_n': np.array([2.5, 3.3, 2.5, 2.5, 2.8, 1.3, 1.5, 1.5, 1.5, 1.5]),
            'eps_n': np.array([1.5, 1.75, 1., 3.75, 8.5, 11.5, 2.5, 3., 2.5, 3.]),

            # Parameters for backscattered (elastically scattered) electrons
            # (25) in FP paper
            'p1EInf': 0.02,      # Minimum probability of elastic scattering (at infinite energy)
            'p1Ehat': 0.496,     # Peak probability
            'eEHat': 0.,        # Peak energy
            'w': 60.86,     # Exponential factor 1
            'p': 1.,        # Exponential factor 2
            # (47a)                 # Angular factors
            'e1': 0.26,
            'e2': 2.,
            # (26)
            'sigmaE': 2.,

            # Parameters for rediffused electrons
            # (28)
            'p1RInf': 0.2,       # Minimum probability of rediffused scattering (at infinite energy)
            'eR': 0.041,     # Peak energy
            'r': 0.104,     # Exponential factor
            # (29)
            'q': 0.5,
            # (47b)                 # Angular factors
            'r1': 0.26,
            'r2': 2.,

            # Parameters for true secondaries
            # (31)
            'deltaTSHat': 1.8848,    # Maximum probability of secondaries
            'eHat0': 276.8,     # Peak enery
            # (32)
            's': 1.54,      # Form factor of fitting curve
            # (48a)                 # Angular factors
            't1': 0.66,
            't2': 0.8,
            # (48b)
            't3': 0.7,
            't4': 1.,
        },
    },
    'combined_simulations_secondaryEmission_machine_parameters': {
        'mandatory': set(),
        'optional': {},
    },
    'additional_cloud_parameters': {
        'mandatory': {

            # Cloud particles
            'cloud_mass',
            'cloud_charge',

            # Residual gas ionization flag
            'gas_ion_flag',

            # Photoemission flag
            'photoem_flag',

            # Uniform initial distribution flag
            'init_unif_flag',

            # Uniform initial density flag
            'init_unif_edens_flag',

            # Secondary emission model
            'switch_model',
        },
        'optional': {

            # MP management settings
            'N_mp_max': (),
            'N_mp_regen': (),
            'N_mp_regen_low': (),
            't_ON_regen_low': (),
            'N_mp_after_regen': (),
            'fact_split': (),
            'fact_clean': (),
            'nel_mp_ref_0': (),
            'Nx_regen': (), 'Ny_regen': (), 'Nvx_regen': (), 'Nvy_regen': (), 'Nvz_regen': (),
            'regen_hist_cut': (),

            'N_mp_soft_regen': (),
            'N_mp_after_soft_regen': (),

            # Tracking and magnetic field
            'N_sub_steps': (),

            # Residual gas ionization
            'P_nTorr': (),
            'sigma_ion_MBarn': (),
            'Temp_K': (),
            'unif_frac': (),
            'E_init_ion': (),

            't_ion': (),

            # Photoemission
            'inv_CDF_refl_photoem_file': (),
            'inv_CDF_all_photoem_file': (),
            'k_pe_st': (),
            'refl_frac': (),
            'alimit': (),
            'e_pe_sigma': (),
            'e_pe_max': (),
            'x0_refl': (),
            'y0_refl': (),
            'out_radius': (),
            'phem_resc_fac': (),
            'photoelectron_angle_distribution': (),
            'energy_distribution': (),
            'flag_continuous_emission': (),
            'filename_chm_photoem': (),

            # Uniform initial distribution
            'Nel_init_unif': (),
            'E_init_unif': (),
            'x_max_init_unif': (),
            'x_min_init_unif': (),
            'y_max_init_unif': (),
            'y_min_init_unif': (),
            'filename_init_MP_state': (),

            # Uniform initial density
            'init_unif_edens': (),
            'E_init_unif_edens': (),
            'x_max_init_unif_edens': (),
            'x_min_init_unif_edens': (),
            'y_max_init_unif_edens': (),
            'y_min_init_unif_edens': (),

            # Secondary emission parameters
            'E_th': (),
            'sigmafit': (),
            'mufit': (),
            'Emax': (),
            's_param': (),
            'del_max': (),
            'R0': (),
            'E0': (),

            'switch_no_increase_energy': (),
            'thresh_low_energy': (),
            'scrub_en_th': (),

            'secondary_angle_distribution': (),

            'sey_file': (),
            'flag_costheta_Emax_shift': (),
            'flag_costheta_delta_scale': (),

            # Saving settings
            'Dx_hist': (),
            'r_center': (),
            'Dt_En_hist': (),
            'Nbin_En_hist': (),
            'En_hist_max': (),

            'flag_movie': (),
            'flag_sc_movie': (),
            'flag_cos_angle_hist': (),
            'cos_angle_width': (),
            'save_mp_state_time_file': (),
            'flag_detailed_MP_info': (),
            'flag_hist_impact_seg': (),
            'flag_verbose_file': (),
            'flag_verbose_stdout': (),
            'dec_fac_secbeam_prof': (),
            'el_density_probes': (),
            'save_simulation_state_time_file': (),
            'x_min_hist_det': (),
            'y_min_hist_det': (),
            'x_max_hist_det': (),
            'y_max_hist_det': (),
            'Dx_hist_det': (),

            # Log and progress files
            'logfile_path': (),
            'progress_path': (),
        },
    },
}

for key in ('secondary_emission_parameters', 'machine_parameters', 'simulation_parameters'):
    parameters_dict['combined_simulations_secondaryEmission_machine_parameters']['mandatory'].update(parameters_dict[key]['mandatory'])
    parameters_dict['combined_simulations_secondaryEmission_machine_parameters']['optional'].update(parameters_dict[key]['optional'])

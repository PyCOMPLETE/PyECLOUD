
parameters_dict = {
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
            'N_mp_soft_regen',

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
        'optional':{
            # Secondary beams
            'secondary_beams_file_list': [],

            'N_mp_soft_regen': None,
            'N_mp_after_soft_regen': None,
            'stopfile': 'stop',

            # Saving settings
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
            'el_density_probes': [],
            'save_simulation_state_time_file': -1,
            'x_min_hist_det': None,
            'y_min_hist_det': None,
            'x_max_hist_det': None,
            'y_max_hist_det': None,
            'Dx_hist_det': None,

            # Undocumented
            'sparse_solver': 'scipy_slu',
            'PyPICmode'    : 'FiniteDifferences_ShortleyWeller',

            # Multigrid parameters
            'f_telescope': None,
            'target_grid': None,
            'N_nodes_discard': None,
            'N_min_Dh_main': None,

            'dec_fact_out': 1,

            # Where to put this?
            't_ion': -1
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

            # Tracking and magnetic field
            'track_method': 'StrongBdip',
            'B': 0., #Tesla (if B=-1 computed from energy and bending radius)
            'bm_totlen': -1, #m
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
            'k_pe_st'                           : -1,
            'refl_frac'                         : -1,
            'alimit'                            : -1,
            'e_pe_sigma'                        : -1,
            'e_pe_max'                          : -1,
            'x0_refl'                           : -1,
            'y0_refl'                           : -1,
            'out_radius'                        : -1,
            'phem_resc_fac'                     : 0.9999,
            'photoelectron_angle_distribution'  : 'cosine_2D',

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
            'E_init_unif_edens'     : None,
            'x_max_init_unif_edens' : None,
            'x_min_init_unif_edens' : None,
            'y_max_init_unif_edens' : None,
            'y_min_init_unif_edens' : None,

            # Undocumented
            'flag_assume_convex': True,
            't_ion': -1 # not in documentation, default value arbitrary
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
            'coast_dens',
            'flag_bunched_beam',
        },
        'optional': {
            # Basic definitions
            'm0_part',
            'Dp_p',
            'nemittx', 'nemitty',
            'x_beam_pos', 'y_beam_pos',
            'sigmax', 'sigmay',

            #Transverse electric field
            'save_beam_field_file_as',

            # if beam_field_file is given
            'Dh_beam_field',
            'Nx', 'Ny',
            'nimag',

            # if compute_FDSW_multigrid
            'Dh_beam_field',
            'f_telescope_beam',
            'target_grid_beam',
            'N_nodes_discard_beam',
            'N_min_Dh_main_beam',

            # if flag_bunched_beam == 1
            'sigmaz',
            't_offs',
            'filling_pattern_file',

            # if flag_bunched_beam == 0
            'beam_long_prof_file',
        },
    },
    'secondary_emission_parameters': {
        'mandatory': {

            # Secondray Electron Yield
            'Emax',
            'del_max',
            'R0',
            'E_th',
            'sigmafit',
            'mufit',

            # Other parameters
            'switch_no_increase_energy',
            'thresh_low_energy',
            'scrub_en_th',
        },
        'optional': {

            # Choice of model
            'switch_model': 0,

            # Undocumented
            'E0': None,
            's_param': None,

            # Other parameters
            'secondary_angle_distribution': 'cosine_2D',
        },
    },
}


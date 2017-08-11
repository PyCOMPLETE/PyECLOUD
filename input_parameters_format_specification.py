import types

class PyECLOUD_ConfigException(Exception):
    pass

def assert_module_has_parameters(module):
    module_name = module.__name__
    essential_parameters = parameters_dict[module_name]['essential']
    optional_parameters = parameters_dict[module_name]['optional']

    module_contents = filter(lambda x: (not x.startswith('__')), dir(module))
    module_contents = set(filter(lambda x: (not isinstance(getattr(module, x), types.ModuleType)), module_contents))

    allowed_parameters = essential_parameters.union(optional_parameters)

    missing_parameters = essential_parameters.difference(module_contents)
    if missing_parameters:
        raise PyECLOUD_ConfigException('Error! These essential parameters are not provided by %s: %r' % (module_name, missing_parameters))

    extra_parameters = module_contents.difference(allowed_parameters)
    if extra_parameters:
        raise PyECLOUD_ConfigException('Error! These parameters should not be in %s: %r' % (module_name, extra_parameters))



parameters_dict = {
    'simulation_parameters': {
        'essential': {

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
            'N_mp_after_soft_regen',

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

            'dec_fact_out',
        },
        'optional':{
            # Secondary beams
            'secondary_beam_file_list',

            'N_mp_soft_regen',
            'stopfile',
            'save_mp_state_time_file',

            # Saving settings
            'flag_movie',
            'flag_sc_movie',
            'save_mp_state_time_file',
            'flag_detailed_MP_info',
            'flag_hist_impact_seg'
            'flag_verbose_file',
            'flag_verbose_stdout',
            'dec_fac_secbeam_prof',
            'el_density_probes',
            'save_simulation_state_time_file',
            'x_min_hist_det', 'y_min_hist_det', 'x_max_hist_det', 'y_max_hist_det', 'Dx_hist_det',

            # Undocumented
            'sparse_solver',

            # Multigrid parameters
            'f_telescope',
            'target_grid',
            'N_nodes_discard',
            'N_min_Dh_main',

            't_ion', # not in documentation
        },
    },
    'machine_parameters': {
        'essential': set(),
        'optional': {

            # Chamber profile
            'chamb_type',
            'x_aper', 'y_aper',
            'filename_chm',

            # Tracking and magnetic field
            'track_method',
            'B',
            'bm_totlen',
            'B_map_file',
            'fact_Bmap',
            'B0x', 'B0y',
            'B_zero_thrhld',
            'N_sub_steps',
            'B_multip',
            'B_skew',

            # Optics
            'betafx', 'betafy',
            'Dx', 'Dy',

            # Residual gas ionization
            'gas_ion_flag',
            'P_nTorr',
            'sigma_ion_MBarn',
            'Temp_K',
            'unif_frac',
            'E_init_ion',
            't_ion', # not in documentation

            # Photoemission
            'photoem_flag',
            'inv_CDF_refl_photoem_file',
            'k_pe_st',
            'refl_frac',
            'alimit',
            'energy_distribution',
            'e_pe_sigma',
            'e_pe_max',
            'x0_refl', 'y0_refl'
            'out_radius',
            'phem_resc_fac',
            'photoelectron_angle_distribution',

            # Uniform initial distribution
            'init_unif_flag',
            'Nel_init_unif',
            'E_init_unif',
            'x_max_init_unif', 'x_min_init_unif', 'y_max_init_unif', 'y_min_init_unif',

        },
    },
    'beam.beam': {
        'essential': {

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
        'essential': set(),
        'optional': {

            # Choice of model
            'switch_model',

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
            'secondary_angle_distribution',
        },
    },
}


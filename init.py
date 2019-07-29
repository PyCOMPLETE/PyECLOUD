#-Begin-preamble-------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 8.1.0
#
#
#     Main author:          Giovanni IADAROLA
#                           BE-ABP Group
#                           CERN
#                           CH-1211 GENEVA 23
#                           SWITZERLAND
#                           giovanni.iadarola@cern.ch
#
#     Contributors:         Eleonora Belli
#                           Philipp Dijkstal
#                           Lotta Mether
#                           Annalisa Romano
#                           Giovanni Rumolo
#                           Eric Wulff
#
#
#     Copyright  CERN,  Geneva  2011  -  Copyright  and  any   other
#     appropriate  legal  protection  of  this  computer program and
#     associated documentation reserved  in  all  countries  of  the
#     world.
#
#     Organizations collaborating with CERN may receive this program
#     and documentation freely and without charge.
#
#     CERN undertakes no obligation  for  the  maintenance  of  this
#     program,  nor responsibility for its correctness,  and accepts
#     no liability whatsoever resulting from its use.
#
#     Program  and documentation are provided solely for the use  of
#     the organization to which they are distributed.
#
#     This program  may  not  be  copied  or  otherwise  distributed
#     without  permission. This message must be retained on this and
#     any other authorized copies.
#
#     The material cannot be sold. CERN should be  given  credit  in
#     all references.
#
#-End-preamble---------------------------------------------------------

from __future__ import division, print_function
import os
import numpy as np
from scipy.constants import c, m_e, e as qe

import myloadmat_to_obj as mlm

import beam_and_timing as beatim
from geom_impact_ellip import ellip_cham_geom_object
import geom_impact_poly_fast_impact as gipfi
import geom_impact_rect_fast_impact as girfi

from sec_emission_model_ECLOUD import SEY_model_ECLOUD
from sec_emission_model_accurate_low_ene import SEY_model_acc_low_ene
from sec_emission_model_ECLOUD_nunif import SEY_model_ECLOUD_non_unif
from sec_emission_model_ECLOUD_nunif import SEY_model_ECLOUD_non_unif_charging
from sec_emission_model_cos_low_ener import SEY_model_cos_le
from sec_emission_model_flat_low_ener import SEY_model_flat_le
from sec_emission_model_from_file import SEY_model_from_file
from sec_emission_model_furman_pivi import SEY_model_furman_pivi
from sec_emission_model_perfect_absorber import SEY_model_perfect_absorber

import dynamics_dipole as dyndip
import dynamics_Boris_f2py as dynB
import dynamics_strong_B_generalized as dyngen
import dynamics_Boris_multipole as dynmul

import MP_system as MPs
import space_charge_class as scc
import impact_management_class as imc
import pyecloud_saver as pysav
import gas_ionization_class as gic
import gen_photoemission_class as gpc

import parse_beam_file as pbf
import parse_cloud_file as pcf
import input_parameters_format_specification as inp_spec
import cloud_manager as cman
import cross_ionization as cion


def read_parameter_files(pyecl_input_folder='./', skip_beam_files=False):
    simulation_param_file = 'simulation_parameters.input'

    # Parse simulation_parameters.input
    input_parameters = inp_spec.import_module_from_file('simulation_parameters', os.path.join(pyecl_input_folder, simulation_param_file))

    # Get names of other input files
    machine_param_file = input_parameters.machine_param_file
    secondary_emission_parameters_file = input_parameters.secondary_emission_parameters_file
    beam_parameters_file = input_parameters.beam_parameters_file

    # Parse other input files
    machine_parameters = inp_spec.import_module_from_file('machine_parameters', os.path.join(pyecl_input_folder, machine_param_file))
    secondary_emission_parameters = inp_spec.import_module_from_file('secondary_emission_parameters', os.path.join(pyecl_input_folder, secondary_emission_parameters_file))

    # Update input_parameters object with parameters from other files
    inp_spec.update_module(input_parameters, machine_parameters)
    inp_spec.update_module(input_parameters, secondary_emission_parameters)
    # Check validity of input files
    inp_spec.assert_module_has_parameters(input_parameters, 'combined_simulations_secondaryEmission_machine_parameters')

    # Create config_dict with all allowed input parameters (not specified are set to default)
    config_dict = {}
    inp_spec.update_config_dict(config_dict, input_parameters, 'combined_simulations_secondaryEmission_machine_parameters')

    # Check validity of main beam input file (not yet used at this stage)
    if not skip_beam_files:
        beam_beam = inp_spec.import_module_from_file('beam_beam', os.path.join(pyecl_input_folder, beam_parameters_file))
        inp_spec.assert_module_has_parameters(beam_beam, 'beam_beam')

    return config_dict


def read_input_files_and_init_components(pyecl_input_folder='./', skip_beam=False,
                                         skip_pyeclsaver=False, skip_spacech_ele=False,
                                         spacech_ele=None,
                                         ignore_kwargs=(), **kwargs):

    config_dict = read_parameter_files(pyecl_input_folder, skip_beam_files=skip_beam)
    # Override config values with kwargs
    for attr, value in kwargs.items():
        if attr in ignore_kwargs:
            continue
        print('Ecloud init. From kwargs: %s = %r' % (attr, value))
        if attr in config_dict:
            config_dict[attr] = value
        else:
            raise inp_spec.PyECLOUD_ConfigException('What exactly does %s do? It is not an expected input.' % attr)

    # Config object
    cc = mlm.obj_from_dict(config_dict)

    # Init beam and possibly second beams
    if not skip_beam:
        flag_presence_sec_beams = len(cc.secondary_beams_file_list) > 0
        b_par = pbf.beam_descr_from_fil(os.path.join(pyecl_input_folder, cc.beam_parameters_file), cc.betafx, cc.Dx, cc.betafy, cc.Dy)

        sec_b_par_list = []
        if flag_presence_sec_beams:
            for sec_b_file in cc.secondary_beams_file_list:
                sec_b_par_list.append(pbf.beam_descr_from_fil(os.path.join(pyecl_input_folder, sec_b_file), cc.betafx, cc.Dx, cc.betafy, cc.Dy))
    else:
        flag_presence_sec_beams = False
        sec_b_par_list = []

    # Init of saver (first print to stdout)
    if not skip_pyeclsaver:
        pyeclsaver = pysav.pyecloud_saver(cc.logfile_path)
    else:
        pyeclsaver = None

    # Parse additional cloud files
    flag_multiple_clouds = len(cc.additional_clouds_file_list) > 0
    print('Simulation with multiple clouds!')

    cloud_par_list = []
    # Make parameter object for default cloud from config dict
    cloud_par_list.append(pcf.cloud_descr_from_file(cloudfilename=None, default_param_obj=cc))
    if flag_multiple_clouds:
        # Get parameter object for additional clouds (unspecified optional parameters set from general config dict)
        for add_cloud_file in cc.additional_clouds_file_list:
            cloud_par_list.append(pcf.cloud_descr_from_file(cloudfilename=os.path.join(pyecl_input_folder, add_cloud_file), default_param_obj=cc))

    # Init chamber
    flag_non_unif_sey = False
    for cloud_par in cloud_par_list:
        if cloud_par.cc.switch_model=="ECLOUD_nunif" or cloud_par.cc.switch_model=="ECLOUD_nunif_charging":
            flag_non_unif_sey = True

    chamber_kwargs = {
        'flag_verbose_file': cc.flag_verbose_file,
        'flag_verbose_stdout': cc.flag_verbose_stdout,
        'flag_assume_convex': cc.flag_assume_convex,
    }

    if cc.chamb_type == 'ellip':
        chamb = ellip_cham_geom_object(cc.x_aper, cc.y_aper, flag_verbose_file=cc.flag_verbose_file)
    elif cc.chamb_type in ('polyg', 'polyg_cython'):
        if os.path.isfile(pyecl_input_folder + '/' + cc.filename_chm):
            filename_chm_path = pyecl_input_folder + '/' + cc.filename_chm
        elif os.path.isfile(pyecl_input_folder + '/' + cc.filename_chm + '.mat'):
            filename_chm_path = pyecl_input_folder + '/' + cc.filename_chm + '.mat'
        else:
            filename_chm_path = cc.filename_chm
        chamb = gipfi.polyg_cham_geom_object(filename_chm_path, flag_non_unif_sey, **chamber_kwargs)
    elif cc.chamb_type == 'rect':
        chamb = girfi.rect_cham_geom_object(cc.x_aper, cc.y_aper, flag_non_unif_sey, **chamber_kwargs)
    else:
        raise inp_spec.PyECLOUD_ConfigException('Chamber type not recognized (choose: ellip/rect/polyg)')

    # Init beam and timing
    if not skip_beam:

        try:
            if os.path.isfile(pyecl_input_folder + '/' + b_par.beam_long_prof_file):
                beam_long_prof_file_path = pyecl_input_folder + '/' + b_par.beam_long_prof_file
            elif os.path.isfile(pyecl_input_folder + '/' + b_par.beam_long_prof_file + '.mat'):
                beam_long_prof_file_path = pyecl_input_folder + '/' + b_par.beam_long_prof_file + '.mat'
            else:
                beam_long_prof_file_path = b_par.beam_long_prof_file
        except:
            beam_long_prof_file_path = b_par.beam_long_prof_file

        if cc.progress_path is not None:
            progress_mapgen_file = cc.progress_path + '_mapgen'
        else:
            progress_mapgen_file = None

        beamtim = beatim.beam_and_timing(b_par.flag_bunched_beam, b_par.fact_beam, b_par.coast_dens, b_par.q_part, b_par.beam_field_file, cc.lam_th,
                                         b_spac=b_par.b_spac, sigmaz=b_par.sigmaz, t_offs=b_par.t_offs, filling_pattern_file=b_par.filling_pattern_file, Dt=cc.Dt, t_end=cc.t_end,
                                         beam_long_prof_file=beam_long_prof_file_path, Dh_beam_field=b_par.Dh_beam_field, f_telescope_beam=b_par.f_telescope_beam,
                                         target_grid_beam=b_par.target_grid_beam, N_nodes_discard_beam=b_par.N_nodes_discard_beam, N_min_Dh_main_beam=b_par.N_min_Dh_main_beam,
                                         chamb=chamb, sigmax=b_par.sigmax, sigmay=b_par.sigmay,
                                         x_beam_pos=b_par.x_beam_pos, y_beam_pos=b_par.y_beam_pos, save_beam_field_file_as=b_par.save_beam_field_file_as,
                                         Nx=b_par.Nx, Ny=b_par.Ny, nimag=b_par.nimag, progress_mapgen_file=progress_mapgen_file)

        sec_beams_list = []
        if flag_presence_sec_beams:
            N_sec_beams = len(sec_b_par_list)
            for ii in xrange(N_sec_beams):
                print('Initialize secondary beam %d/%d' % (ii + 1, N_sec_beams))
                sb_par = sec_b_par_list[ii]

                try:
                    if os.path.isfile(pyecl_input_folder + '/' + sb_par.beam_long_prof_file):
                        sbeam_long_prof_file_path = pyecl_input_folder + '/' + sb_par.beam_long_prof_file
                    elif os.path.isfile(pyecl_input_folder + '/' + b_par.beam_long_prof_file + '.mat'):
                        sbeam_long_prof_file_path = pyecl_input_folder + '/' + sb_par.beam_long_prof_file + '.mat'
                    else:
                        sbeam_long_prof_file_path = sb_par.beam_long_prof_file
                except TypeError:
                    # in case sb_par.beam_long_prof_file is -1
                    sbeam_long_prof_file_path = sb_par.beam_long_prof_file

                sec_beams_list.append(beatim.beam_and_timing(sb_par.flag_bunched_beam, sb_par.fact_beam, sb_par.coast_dens, sb_par.q_part, sb_par.beam_field_file, cc.lam_th,
                                                             b_spac=sb_par.b_spac, sigmaz=sb_par.sigmaz, t_offs=sb_par.t_offs, filling_pattern_file=sb_par.filling_pattern_file, Dt=cc.Dt, t_end=cc.t_end,
                                                             beam_long_prof_file=sbeam_long_prof_file_path, Dh_beam_field=sb_par.Dh_beam_field, f_telescope_beam=sb_par.f_telescope_beam,
                                                             target_grid_beam=sb_par.target_grid_beam, N_nodes_discard_beam=sb_par.N_nodes_discard_beam, N_min_Dh_main_beam=sb_par.N_min_Dh_main_beam,
                                                             chamb=chamb, sigmax=sb_par.sigmax, sigmay=sb_par.sigmay,
                                                             x_beam_pos=sb_par.x_beam_pos, y_beam_pos=sb_par.y_beam_pos, save_beam_field_file_as=sb_par.save_beam_field_file_as,
                                                             flag_secodary_beam=True, t_primary_beam=beamtim.t,
                                                             Nx=sb_par.Nx, Ny=sb_par.Ny, nimag=sb_par.nimag, progress_mapgen_file=(cc.progress_path + ('_mapgen_sec_%d' % ii))))
    else:
        beamtim = None
        sec_beams_list = []

    # Init spacecharge
    if skip_spacech_ele:
        spacech_ele_sim = None
    elif spacech_ele is not None:
        spacech_ele_sim = spacech_ele
    else:
        if cc.sparse_solver == 'klu':
            print('''sparse_solver: 'klu' no longer supported --> going to PyKLU''')
            cc.sparse_solver = 'PyKLU'
        spacech_ele_sim = scc.space_charge(chamb, cc.Dh_sc, Dt_sc=cc.Dt_sc, sparse_solver=cc.sparse_solver, PyPICmode=cc.PyPICmode,
                                           f_telescope=cc.f_telescope, target_grid=cc.target_grid, N_nodes_discard=cc.N_nodes_discard, N_min_Dh_main=cc.N_min_Dh_main)

    flag_cross_ion = False
    if cc.cross_ion_definitions is not None:
        flag_cross_ion = True

    # Loop over clouds to init all cloud-specific objects
    cloud_list = []
    for cloud_par in cloud_par_list:
        thiscloud = cloud_par.cc

        print('Initialize cloud %s:' % (thiscloud.cloud_name))

        # Init saver for all but default cloud (which is already initialized)
        if cloud_par is not cloud_par_list[0]:
            if not skip_pyeclsaver:
                pyeclsaver = pysav.pyecloud_saver(thiscloud.logfile_path)
            else:
                pyeclsaver = None

        # Init MP system
        MP_e = MPs.MP_system(thiscloud.N_mp_max, thiscloud.nel_mp_ref_0, thiscloud.fact_split, thiscloud.fact_clean,
                             thiscloud.N_mp_regen_low, thiscloud.N_mp_regen, thiscloud.N_mp_after_regen,
                             thiscloud.Dx_hist, thiscloud.Nx_regen, thiscloud.Ny_regen, thiscloud.Nvx_regen,
                             thiscloud.Nvy_regen, thiscloud.Nvz_regen, thiscloud.regen_hist_cut, chamb,
                             N_mp_soft_regen=thiscloud.N_mp_soft_regen, N_mp_after_soft_regen=thiscloud.N_mp_after_soft_regen,
                             N_mp_async_regen=thiscloud.N_mp_async_regen, N_mp_after_async_regen=thiscloud.N_mp_after_async_regen,
                             charge=thiscloud.cloud_charge, mass=thiscloud.cloud_mass, flag_lifetime_hist = thiscloud.flag_lifetime_hist,
                             name=thiscloud.cloud_name)

        # Init secondary emission object
        if thiscloud.switch_model == 'perfect_absorber':
            sey_mod = SEY_model_perfect_absorber()
        else:

            kwargs_secem = {}
            if thiscloud.E0 is not None:
                kwargs_secem.update({'E0': thiscloud.E0})
                #If E0 is not provided use default value for each object
            if thiscloud.s_param is not None:
                if thiscloud.switch_model == 0 or thiscloud.switch_model == 'ECLOUD':
                    kwargs_secem.update({'s': thiscloud.s_param})
                else:
                    raise inp_spec.PyECLOUD_ConfigException('s parameter can be changed only in the ECLOUD sec. emission model!')

            if thiscloud.switch_model in (0, 'ECLOUD'):
                kwargs_secem['flag_costheta_delta_scale'] = thiscloud.flag_costheta_delta_scale
                kwargs_secem['flag_costheta_Emax_shift'] = thiscloud.flag_costheta_Emax_shift
                sey_mod = SEY_model_ECLOUD(
                    thiscloud.Emax, thiscloud.del_max, thiscloud.R0,
                    E_th=thiscloud.E_th, sigmafit=thiscloud.sigmafit, mufit=thiscloud.mufit,
                    switch_no_increase_energy=thiscloud.switch_no_increase_energy,
                    thresh_low_energy=thiscloud.thresh_low_energy,
                    secondary_angle_distribution=thiscloud.secondary_angle_distribution,
                    **kwargs_secem)
            elif thiscloud.switch_model in (1, 'ACC_LOW'):
                sey_mod = SEY_model_acc_low_ene(thiscloud.Emax, thiscloud.del_max, thiscloud.R0,
                                                E_th=thiscloud.E_th, sigmafit=thiscloud.sigmafit, mufit=thiscloud.mufit,
                                                switch_no_increase_energy=thiscloud.switch_no_increase_energy,
                                                thresh_low_energy=thiscloud.thresh_low_energy,
                                                secondary_angle_distribution=thiscloud.secondary_angle_distribution,
                                                **kwargs_secem)
            elif thiscloud.switch_model == 'ECLOUD_nunif':
                sey_mod = SEY_model_ECLOUD_non_unif(chamb, thiscloud.Emax, thiscloud.del_max, thiscloud.R0,
                                                    E_th=thiscloud.E_th, sigmafit=thiscloud.sigmafit, mufit=thiscloud.mufit,
                                                    switch_no_increase_energy=thiscloud.switch_no_increase_energy,
                                                    thresh_low_energy=thiscloud.thresh_low_energy,
                                                    secondary_angle_distribution=thiscloud.secondary_angle_distribution,
                                                    **kwargs_secem)
            elif thiscloud.switch_model == 'ECLOUD_nunif_charging':
                sey_mod = SEY_model_ECLOUD_non_unif_charging(chamb, thiscloud.Emax, thiscloud.del_max, thiscloud.R0,
                                                    E_th=thiscloud.E_th, sigmafit=thiscloud.sigmafit, mufit=thiscloud.mufit,
                                                    switch_no_increase_energy=thiscloud.switch_no_increase_energy,
                                                    thresh_low_energy=thiscloud.thresh_low_energy,
                                                    secondary_angle_distribution=thiscloud.secondary_angle_distribution,
                                                    **kwargs_secem)
            elif thiscloud.switch_model == 'cos_low_ene':
                sey_mod = SEY_model_cos_le(thiscloud.Emax, thiscloud.del_max, thiscloud.R0,
                                           E_th=thiscloud.E_th, sigmafit=thiscloud.sigmafit, mufit=thiscloud.mufit,
                                           switch_no_increase_energy=thiscloud.switch_no_increase_energy,
                                           thresh_low_energy=thiscloud.thresh_low_energy,
                                           secondary_angle_distribution=thiscloud.secondary_angle_distribution,
                                           **kwargs_secem)
            elif thiscloud.switch_model == 'flat_low_ene':
                sey_mod = SEY_model_flat_le(thiscloud.Emax, thiscloud.del_max, thiscloud.R0,
                                            E_th=thiscloud.E_th, sigmafit=thiscloud.sigmafit, mufit=thiscloud.mufit,
                                            switch_no_increase_energy=thiscloud.switch_no_increase_energy,
                                            thresh_low_energy=thiscloud.thresh_low_energy,
                                            secondary_angle_distribution=thiscloud.secondary_angle_distribution,
                                            **kwargs_secem)
            elif thiscloud.switch_model == 'from_file':
                kwargs_secem['flag_costheta_delta_scale'] = thiscloud.flag_costheta_delta_scale
                kwargs_secem['flag_costheta_Emax_shift'] = thiscloud.flag_costheta_Emax_shift
                if os.path.isfile(pyecl_input_folder + '/' + thiscloud.sey_file):
                    sey_file_path = pyecl_input_folder + '/' + thiscloud.sey_file
                elif os.path.isfile(pyecl_input_folder + '/' + thiscloud.sey_file + '.mat'):
                    sey_file_path = pyecl_input_folder + '/' + thiscloud.sey_file + '.mat'
                else:
                    sey_file_path = thiscloud.sey_file
                sey_mod = SEY_model_from_file(sey_file_path,
                                              E_th=thiscloud.E_th, sigmafit=thiscloud.sigmafit, mufit=thiscloud.mufit,
                                              switch_no_increase_energy=thiscloud.switch_no_increase_energy,
                                              thresh_low_energy=thiscloud.thresh_low_energy,
                                              secondary_angle_distribution=thiscloud.secondary_angle_distribution,
                                              **kwargs_secem)
            elif(thiscloud.switch_model == 'furman_pivi'):
                kwargs_secem['flag_costheta_delta_scale'] = thiscloud.flag_costheta_delta_scale
                kwargs_secem['flag_costheta_Emax_shift'] = thiscloud.flag_costheta_Emax_shift
                sey_mod = SEY_model_furman_pivi(E_th=thiscloud.E_th, sigmafit=thiscloud.sigmafit, mufit=thiscloud.mufit,
                                                switch_no_increase_energy=thiscloud.switch_no_increase_energy,
                                                thresh_low_energy=thiscloud.thresh_low_energy,
                                                secondary_angle_distribution=thiscloud.secondary_angle_distribution,
                                                furman_pivi_surface=thiscloud.furman_pivi_surface,
                                                **kwargs_secem)

            else:
                raise inp_spec.PyECLOUD_ConfigException('switch_model not recognized!')

        # Init impact management
        flag_seg = (thiscloud.flag_hist_impact_seg == 1 or thiscloud.flag_hist_impact_seg is True)

        if flag_seg and cc.chamb_type == 'ellip':
            print('Warning: You cannot enable flag_hist_impact_seg for an ellip chamber --> disabled!')
            flag_seg = False

        if cc.flag_lifetime_hist:
            if cc.Nbin_lifetime_hist is None or cc.lifetime_hist_max is None or cc.Dt_lifetime_hist is None:
                raise inp_spec.PyECLOUD_ConfigException(
                        'If flag_lifetime_hist is True then all the histogram parameters must be specified')

        impact_man = imc.impact_management(chamb, sey_mod,
            thiscloud.Dx_hist, thiscloud.scrub_en_th, cc.Nbin_En_hist, cc.En_hist_max,
            cc.Nbin_lifetime_hist, cc.lifetime_hist_max, cc.flag_lifetime_hist,
            flag_seg=flag_seg, flag_En_hist_seg=thiscloud.flag_En_hist_seg,
            cos_angle_width=cc.cos_angle_width)

        # Init gas ionization and photoemission
        if thiscloud.gas_ion_flag == 1:
            resgasion = gic.residual_gas_ionization(thiscloud.unif_frac, thiscloud.P_nTorr, thiscloud.sigma_ion_MBarn,
                                                    thiscloud.Temp_K, chamb, thiscloud.E_init_ion, thiscloud.flag_lifetime_hist)
        else:
            resgasion = None

        if thiscloud.photoem_flag == 1:
            phemiss = gpc.photoemission(thiscloud.inv_CDF_refl_photoem_file, thiscloud.k_pe_st, thiscloud.refl_frac, thiscloud.e_pe_sigma, thiscloud.e_pe_max,
                                        thiscloud.alimit, thiscloud.x0_refl, thiscloud.y0_refl, thiscloud.out_radius, chamb, thiscloud.phem_resc_fac,
                                        thiscloud.energy_distribution, thiscloud.photoelectron_angle_distribution, beamtim, thiscloud.flag_continuous_emission)
        elif thiscloud.photoem_flag in (2, 'from_file'):
            phemiss = gpc.photoemission_from_file(thiscloud.inv_CDF_all_photoem_file, chamb, thiscloud.phem_resc_fac, thiscloud.energy_distribution,
                                                  thiscloud.e_pe_sigma, thiscloud.e_pe_max, thiscloud.k_pe_st, thiscloud.out_radius,
                                                  thiscloud.photoelectron_angle_distribution, beamtim, thiscloud.flag_continuous_emission)
        elif thiscloud.photoem_flag in (3, 'per_segment'):

            if os.path.isfile(pyecl_input_folder + '/' + thiscloud.filename_chm_photoem):
                filename_chm_photoem_path = pyecl_input_folder + '/' + thiscloud.filename_chm_photoem
            elif os.path.isfile(pyecl_input_folder + '/' + thiscloud.filename_chm_photoem + '.mat'):
                filename_chm_photoem_path = pyecl_input_folder + '/' + thiscloud.filename_chm_photoem + '.mat'
            else:
                filename_chm_photoem_path = thiscloud.filename_chm_photoem

            chamb_phemiss = gipfi.polyg_cham_photoemission(filename_chm_photoem_path)
            if not chamb_phemiss.vertexes_are_subset(chamb):
                raise gipfi.PyECLOUD_ChamberException('Chambers for secondary emission and photoemission do not have the same shape!')
            phemiss = gpc.photoemission_per_segment(chamb_phemiss, thiscloud.energy_distribution, thiscloud.e_pe_sigma, thiscloud.e_pe_max, thiscloud.k_pe_st,
                                                    thiscloud.photoelectron_angle_distribution, beamtim, thiscloud.flag_continuous_emission)
        else:
            phemiss = None

        # Real saver init
        if not skip_pyeclsaver:
            flag_last_cloud = cloud_par is cloud_par_list[-1]
            pyeclsaver.start_observing(cc.Dt, MP_e, beamtim, impact_man,
                                       thiscloud.r_center, thiscloud.Dt_En_hist, thiscloud.logfile_path, thiscloud.progress_path,
                                       flag_detailed_MP_info=thiscloud.flag_detailed_MP_info, flag_movie=thiscloud.flag_movie,
                                       flag_sc_movie=thiscloud.flag_sc_movie, save_mp_state_time_file=thiscloud.save_mp_state_time_file,
                                       flag_presence_sec_beams=flag_presence_sec_beams, sec_beams_list=sec_beams_list, dec_fac_secbeam_prof=thiscloud.dec_fac_secbeam_prof,
                                       el_density_probes=thiscloud.el_density_probes, save_simulation_state_time_file=thiscloud.save_simulation_state_time_file,
                                       x_min_hist_det=thiscloud.x_min_hist_det, x_max_hist_det=thiscloud.x_max_hist_det,
                                       y_min_hist_det=thiscloud.y_min_hist_det, y_max_hist_det=thiscloud.y_max_hist_det,
                                       Dx_hist_det=thiscloud.Dx_hist_det, dec_fact_out=cc.dec_fact_out, stopfile=cc.stopfile, filen_main_outp=thiscloud.filen_main_outp,
                                       flag_cos_angle_hist=thiscloud.flag_cos_angle_hist, cos_angle_width=thiscloud.cos_angle_width,
                                       flag_multiple_clouds=flag_multiple_clouds, cloud_name=thiscloud.cloud_name, flag_last_cloud=flag_last_cloud,
                                       checkpoint_DT=cc.checkpoint_DT, checkpoint_folder=cc.checkpoint_folder, copy_main_outp_folder=cc.copy_main_outp_folder,
                                       copy_main_outp_DT=cc.copy_main_outp_DT, extract_sey=cc.extract_sey,
                                       step_by_step_custom_observables=cc.step_by_step_custom_observables,
                                       pass_by_pass_custom_observables=cc.pass_by_pass_custom_observables,
                                       save_once_custom_observables=cc.save_once_custom_observables, 
                                       flag_lifetime_hist = thiscloud.flag_lifetime_hist,
                                       Dt_lifetime_hist = thiscloud.Dt_lifetime_hist,
                                       extract_ene_dist=cc.extract_ene_dist,
                                       ene_dist_test_E_impact_eV=cc.ene_dist_test_E_impact_eV,
                                       Nbin_extract_ene=cc.Nbin_extract_ene,
                                       factor_ene_dist_max=cc.factor_ene_dist_max,
                                       flag_cross_ion=flag_cross_ion,
                                       save_only = thiscloud.save_only
                                       )
            print('pyeclsaver saves to file: %s' % pyeclsaver.filen_main_outp)

        # Init electron tracker
        if cc.track_method == 'Boris':
            dynamics = dynB.pusher_Boris(cc.Dt, cc.B0x, cc.B0y, cc.B0z,
                                         cc.B_map_file, cc.fact_Bmap, cc.Bz_map_file, N_sub_steps=thiscloud.N_sub_steps)
        elif cc.track_method == 'StrongBdip':
            #~ raise ValueError('The StrongBdip tracker is no longer supported! If you really want to use it remove this line.')
            if not(np.abs(thiscloud.cloud_charge - (-qe)) / np.abs(qe) < 1e-3 and np.abs(thiscloud.cloud_mass - m_e) / m_e < 1e-3):
                raise ValueError('StrongBdip tracking method is implemented only for electrons!')
            if cc.B == -1:
                B = 2 * np.pi * b_par.beta_rel * b_par.energy_J / (c * qe * cc.bm_totlen)
            else:
                B = cc.B
            dynamics = dyndip.pusher_dipole_magnet(cc.Dt, B)
        elif cc.track_method == 'StrongBgen':
            #~ raise ValueError('The StrongBgen tracker is no longer supported! If you really want to use it remove this line.')
            if not(np.abs(thiscloud.cloud_charge - (-qe)) / np.abs(qe) < 1e-3 and np.abs(thiscloud.cloud_mass - m_e) / m_e < 1e-3):
                raise ValueError('StrongBgen tracking method is implemented only for electrons!')
            dynamics = dyngen.pusher_strong_B_generalized(cc.Dt, cc.B0x, cc.B0y,
                                                          cc.B_map_file, cc.fact_Bmap, cc.B_zero_thrhld)
        elif cc.track_method == 'BorisMultipole':
            dynamics = dynmul.pusher_Boris_multipole(Dt=cc.Dt, N_sub_steps=cc.N_sub_steps, B_multip=cc.B_multip, B_skew=cc.B_skew)
        else:
            raise inp_spec.PyECLOUD_ConfigException("track_method should be 'Boris' or 'StrongBdip' or 'StrongBgen' or 'BorisMultipole'")

        # Initial electron density
        if thiscloud.init_unif_flag == 1:
            print("Adding inital %.2e electrons to the initial distribution" % thiscloud.Nel_init_unif)
            MP_e.add_uniform_MP_distrib(thiscloud.Nel_init_unif, thiscloud.E_init_unif,
                                        thiscloud.x_max_init_unif, thiscloud.x_min_init_unif,
                                        thiscloud.y_max_init_unif, thiscloud.y_min_init_unif)

        if thiscloud.init_unif_edens_flag == 1:
            print("Adding inital %.2e electrons/m^3 to the initial distribution" % thiscloud.init_unif_edens)
            MP_e.add_uniform_ele_density(n_ele=thiscloud.init_unif_edens, E_init=thiscloud.E_init_unif_edens,
                                         x_max=thiscloud.x_max_init_unif_edens, x_min=thiscloud.x_min_init_unif_edens,
                                         y_max=thiscloud.y_max_init_unif_edens, y_min=thiscloud.y_min_init_unif_edens)

        if thiscloud.filename_init_MP_state != -1 and thiscloud.filename_init_MP_state is not None:
            print("Adding inital electrons from: %s" % thiscloud.filename_init_MP_state)
            MP_e.add_from_file(thiscloud.filename_init_MP_state)

        # Init empty rho for cloud
        if hasattr(spacech_ele_sim, 'rho'):
            rho = spacech_ele_sim.rho * 0.
        else:
            rho = None

        cloud = cman.Cloud(thiscloud.cloud_name, thiscloud, MP_e, impact_man, dynamics, pyeclsaver, thiscloud.gas_ion_flag,
                           resgasion, thiscloud.t_ion, thiscloud.photoem_flag, phemiss, rho)

        cloud_list.append(cloud)

    # Init cross-ionization
    if flag_cross_ion:
        cross_ion = cion.Cross_Ionization(pyecl_input_folder, cc.cross_ion_definitions, cloud_list)
    else:
        cross_ion = None


    return (beamtim,
            spacech_ele_sim,
            cc.t_sc_ON,
            flag_presence_sec_beams,
            sec_beams_list,
            config_dict,
            flag_multiple_clouds,
            cloud_list,
            cc.checkpoint_folder,
            cross_ion
            )

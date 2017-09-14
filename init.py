#----------------------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                          PyECLOUD Version 6.4.0
#
#
#     Author and contact:   Giovanni IADAROLA
#                           BE-ABP Group
#                           CERN
#                           CH-1211 GENEVA 23
#                           SWITZERLAND
#                           giovanni.iadarola@cern.ch
#
#                contact:   Giovanni RUMOLO
#                           BE-ABP Group
#                           CERN
#                           CH-1211 GENEVA 23
#                           SWITZERLAND
#                           giovanni.rumolo@cern.ch
#
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
#----------------------------------------------------------------------


import numpy as np
from scipy.constants import c, e as qe

from . import myloadmat_to_obj as mlm

from . import beam_and_timing as beatim
from .geom_impact_ellip import ellip_cham_geom_object

from .sec_emission_model_ECLOUD import SEY_model_ECLOUD
from .sec_emission_model_accurate_low_ene import SEY_model_acc_low_ene
from .sec_emission_model_ECLOUD_nunif import SEY_model_ECLOUD_non_unif
from .sec_emission_model_cos_low_ener import SEY_model_cos_le
from .sec_emission_model_flat_low_ener import SEY_model_flat_le

from . import dynamics_dipole as dyndip
from . import dynamics_Boris_f2py as dynB
from . import dynamics_strong_B_generalized as dyngen

from . import MP_system as MPs
from . import space_charge_class as scc
from . import impact_management_class as imc
from . import pyecloud_saver as pysav
from . import gas_ionization_class as gic
from . import gen_photoemission_class as gpc

from . import parse_beam_file as pbf
from . import input_parameters_format_specification as inp_spec


def read_input_files_and_init_components(pyecl_input_folder='./', **kwargs):

    simulation_param_file = 'simulation_parameters.input'
    config_dict = {}

    simulation_parameters = inp_spec.import_module_from_file('simulation_parameters', pyecl_input_folder+'/'+simulation_param_file)
    inp_spec.assert_module_has_parameters(simulation_parameters, 'simulation_parameters')
    inp_spec.update_config_dict(config_dict, simulation_parameters, 'simulation_parameters')

    machine_param_file = config_dict['machine_param_file']
    secondary_emission_parameters_file = config_dict['secondary_emission_parameters_file']
    beam_parameters_file = config_dict['beam_parameters_file']

    machine_parameters = inp_spec.import_module_from_file('machine_parameters', pyecl_input_folder+'/'+machine_param_file)
    inp_spec.assert_module_has_parameters(machine_parameters, 'machine_parameters')
    inp_spec.update_config_dict(config_dict, machine_parameters, 'machine_parameters')

    secondary_emission_parameters = inp_spec.import_module_from_file('secondary_emission_parameters', pyecl_input_folder+'/'+secondary_emission_parameters_file)
    inp_spec.assert_module_has_parameters(secondary_emission_parameters, 'secondary_emission_parameters')
    inp_spec.update_config_dict(config_dict, secondary_emission_parameters, 'secondary_emission_parameters')

    beam_beam = inp_spec.import_module_from_file('beam_beam', pyecl_input_folder+'/'+beam_parameters_file)
    inp_spec.assert_module_has_parameters(beam_beam, 'beam_beam')

    # Probably promote this to an optional parameter in simulation_parameters.input
    filen_main_outp = 'Pyecltest'

    # Override config values with kwargs
    for attr, value in list(kwargs.items()):
        print('Ecloud init. From kwargs: %s = %r' % (attr, value))
        if attr == 'filen_main_outp':
            filen_main_outp = value
        elif attr in config_dict:
            config_dict[attr] = value
        else:
            print('Warning! What exactly does %s do? It is not part of any config file.' % attr)
            exec(('%s=value') % attr)

    cc = mlm.obj_from_dict(config_dict)

    flag_presence_sec_beams = len(cc.secondary_beams_file_list)>0
    b_par = pbf.beam_descr_from_fil(pyecl_input_folder+'/'+cc.beam_parameters_file, cc.betafx, cc.Dx, cc.betafy, cc.Dy)

    sec_b_par_list=[]
    if flag_presence_sec_beams:
        for sec_b_file in cc.secondary_beams_file_list:
            sec_b_par_list.append(pbf.beam_descr_from_fil(pyecl_input_folder+'/'+cc.sec_b_file, cc.betafx, cc.Dx, cc.betafy, cc.Dy))

    if cc.B==-1:
        cc.B = 2*np.pi*b_par.beta_rel*b_par.energy_J/(c*qe*cc.bm_totlen)


    ##########################################

    pyeclsaver=pysav.pyecloud_saver(cc.logfile_path)

    if cc.switch_model=='ECLOUD_nunif':
        flag_non_unif_sey = 1
    else:
        flag_non_unif_sey = 0

    if cc.chamb_type=='ellip':
        chamb=ellip_cham_geom_object(cc.x_aper, cc.y_aper, flag_verbose_file=cc.flag_verbose_file)
    elif cc.chamb_type in ('polyg', 'polyg_cython'):

        from . import geom_impact_poly_fast_impact as gipfi
        chamb=gipfi.polyg_cham_geom_object(cc.filename_chm, flag_non_unif_sey,
                                     flag_verbose_file=cc.flag_verbose_file, flag_verbose_stdout=cc.flag_verbose_stdout, flag_assume_convex=cc.flag_assume_convex)
    elif cc.chamb_type=='polyg_numpy':
        raise ValueError("chamb_type='polyg_numpy' not supported anymore")
        #~ chamb=gip.polyg_cham_geom_object(filename_chm, flag_non_unif_sey,
        #~ flag_verbose_file=flag_verbose_file, flag_verbose_stdout=flag_verbose_stdout)
    elif cc.chamb_type=='rect':
        from . import geom_impact_rect_fast_impact as girfi
        chamb = girfi.rect_cham_geom_object(cc.x_aper, cc.y_aper, flag_verbose_file=cc.flag_verbose_file, flag_verbose_stdout=cc.flag_verbose_stdout)
    else:
        raise ValueError('Chamber type not recognized (choose: ellip/rect/polyg)')


    MP_e=MPs.MP_system(cc.N_mp_max, cc.nel_mp_ref_0, cc.fact_split, cc.fact_clean,
                       cc.N_mp_regen_low, cc.N_mp_regen, cc.N_mp_after_regen,
                       cc.Dx_hist, cc.Nx_regen, cc.Ny_regen, cc.Nvx_regen, cc.Nvy_regen, cc.Nvz_regen, cc.regen_hist_cut, chamb,
                       N_mp_soft_regen=cc.N_mp_soft_regen, N_mp_after_soft_regen=cc.N_mp_after_soft_regen)

    beamtim=beatim.beam_and_timing(b_par.flag_bunched_beam, b_par.fact_beam, b_par.coast_dens, b_par.beam_field_file,cc.lam_th,
                 b_spac=b_par.b_spac, sigmaz=b_par.sigmaz,t_offs=b_par.t_offs, filling_pattern_file=b_par.filling_pattern_file, Dt=cc.Dt, t_end=cc.t_end,
                 beam_long_prof_file=b_par.beam_long_prof_file, Dh_beam_field=b_par.Dh_beam_field, f_telescope_beam=b_par.f_telescope_beam,
                 target_grid_beam=b_par.target_grid_beam, N_nodes_discard_beam=b_par.N_nodes_discard_beam, N_min_Dh_main_beam=b_par.N_min_Dh_main_beam,
                 chamb=chamb,  sigmax=b_par.sigmax, sigmay=b_par.sigmay,
                 x_beam_pos=b_par.x_beam_pos, y_beam_pos=b_par.y_beam_pos, save_beam_field_file_as=b_par.save_beam_field_file_as,
                 Nx=b_par.Nx, Ny=b_par.Ny, nimag=b_par.nimag, progress_mapgen_file=(cc.progress_path+'_mapgen'))

    if cc.sparse_solver=='klu':
        print('''sparse_solver: 'klu' no longer supported --> going to PyKLU''')
        cc.sparse_solver='PyKLU'

    spacech_ele = scc.space_charge(chamb, cc.Dh_sc, Dt_sc=cc.Dt_sc, sparse_solver=cc.sparse_solver, PyPICmode=cc.PyPICmode,
                        f_telescope=cc.f_telescope, target_grid=cc.target_grid, N_nodes_discard=cc.N_nodes_discard, N_min_Dh_main=cc.N_min_Dh_main)

    sec_beams_list=[]
    if flag_presence_sec_beams:
        N_sec_beams = len(sec_b_par_list)
        for ii in range(N_sec_beams):
            print('Initialize secondary beam %d/%d' % (ii+1, N_sec_beams))
            sb_par = sec_b_par_list[ii]
            sec_beams_list.append(beatim.beam_and_timing(sb_par.flag_bunched_beam, sb_par.fact_beam, sb_par.coast_dens, sb_par.beam_field_file, cc.lam_th,
                 b_spac=sb_par.b_spac, sigmaz=sb_par.sigmaz,t_offs=sb_par.t_offs, filling_pattern_file=sb_par.filling_pattern_file, Dt=cc.Dt, t_end=cc.t_end,
                 beam_long_prof_file=sb_par.beam_long_prof_file, Dh_beam_field=sb_par.Dh_beam_field, f_telescope_beam=sb_par.f_telescope_beam,
                 target_grid_beam=sb_par.target_grid_beam, N_nodes_discard_beam=sb_par.N_nodes_discard_beam, N_min_Dh_main_beam=sb_par.N_min_Dh_main_beam,
                 chamb=chamb, sigmax=sb_par.sigmax, sigmay=sb_par.sigmay,
                 x_beam_pos=sb_par.x_beam_pos, y_beam_pos=sb_par.y_beam_pos, save_beam_field_file_as=sb_par.save_beam_field_file_as,
                 flag_secodary_beam=True, t_primary_beam=beamtim.t,
                 Nx=sb_par.Nx, Ny=sb_par.Ny, nimag=sb_par.nimag, progress_mapgen_file=(cc.progress_path+('_mapgen_sec_%d' % ii))))

    kwargs = {}

    if cc.E0 is not None:
        kwargs.update({'E0':cc.E0})
        #If E0 is not provided use default value for each object

    if cc.s_param is not None:
        if cc.switch_model==0 or cc.switch_model=='ECLOUD':
            kwargs.update({'s':cc.s_param})
        else:
            raise ValueError('s parameter can be changed only in the ECLOUD sec. emission model!')

    if cc.switch_model==0 or cc.switch_model=='ECLOUD':
        sey_mod=SEY_model_ECLOUD(cc.Emax,cc.del_max,cc.R0,**kwargs)
    elif cc.switch_model==1 or cc.switch_model=='ACC_LOW':
        sey_mod=SEY_model_acc_low_ene(cc.Emax,cc.del_max,cc.R0,**kwargs)
    elif cc.switch_model=='ECLOUD_nunif':
        sey_mod=SEY_model_ECLOUD_non_unif(chamb, cc.Emax,cc.del_max,cc.R0,**kwargs)
    elif cc.switch_model=='cos_low_ene':
        sey_mod=SEY_model_cos_le(cc.Emax,cc.del_max,cc.R0,**kwargs)
    elif cc.switch_model=='flat_low_ene':
        sey_mod=SEY_model_flat_le(cc.Emax,cc.del_max,cc.R0)


    flag_seg = (cc.flag_hist_impact_seg==1)

    impact_man=imc.impact_management(cc.switch_no_increase_energy, chamb, sey_mod, cc.E_th, cc.sigmafit, cc.mufit,
                 cc.Dx_hist, cc.scrub_en_th, cc.Nbin_En_hist, cc.En_hist_max, thresh_low_energy=cc.thresh_low_energy,
                 flag_seg=flag_seg, cos_angle_width=cc.cos_angle_width, secondary_angle_distribution=cc.secondary_angle_distribution)

    #resgasion_sec_beam_list=[]
    if cc.gas_ion_flag==1:
        resgasion=gic.residual_gas_ionization(cc.unif_frac, cc.P_nTorr, cc.sigma_ion_MBarn,cc.Temp_K,chamb,cc.E_init_ion)
    else:
        resgasion=None


    if cc.photoem_flag:
        phemiss=gpc.photoemission(cc.inv_CDF_refl_photoem_file, cc.k_pe_st, cc.refl_frac, cc.e_pe_sigma, cc.e_pe_max, cc.alimit,
                                  cc.x0_refl, cc.y0_refl, cc.out_radius, chamb, cc.phem_resc_fac, cc.photoelectron_angle_distribution)
    else:
        phemiss=None

    pyeclsaver.start_observing(MP_e, beamtim, impact_man,
                 cc.r_center, cc.Dt_En_hist, cc.logfile_path, cc.progress_path, flag_detailed_MP_info=cc.flag_detailed_MP_info,
                 flag_movie=cc.flag_movie, flag_sc_movie=cc.flag_sc_movie, save_mp_state_time_file=cc.save_mp_state_time_file,
                 flag_presence_sec_beams=flag_presence_sec_beams, sec_beams_list=sec_beams_list, dec_fac_secbeam_prof=cc.dec_fac_secbeam_prof,
                 el_density_probes=cc.el_density_probes, save_simulation_state_time_file=cc.save_simulation_state_time_file,
                 x_min_hist_det=cc.x_min_hist_det, x_max_hist_det=cc.x_max_hist_det, y_min_hist_det=cc.y_min_hist_det, y_max_hist_det=cc.y_max_hist_det,
                 Dx_hist_det=cc.Dx_hist_det, dec_fact_out=cc.dec_fact_out, stopfile=cc.stopfile, filen_main_outp=filen_main_outp,
                 flag_cos_angle_hist=cc.flag_cos_angle_hist, cos_angle_width=cc.cos_angle_width)


    if cc.track_method == 'Boris':
        dynamics=dynB.pusher_Boris(cc.Dt, cc.B0x, cc.B0y, cc.B0z,
                 cc.B_map_file, cc.fact_Bmap, cc.Bz_map_file, N_sub_steps=cc.N_sub_steps)
    elif cc.track_method == 'StrongBdip':
        dynamics=dyndip.pusher_dipole_magnet(cc.Dt,cc.B)
    elif cc.track_method == 'StrongBgen':
        dynamics=dyngen.pusher_strong_B_generalized(cc.Dt, cc.B0x, cc.B0y,
                 cc.B_map_file, cc.fact_Bmap, cc.B_zero_thrhld)
    elif cc.track_method == 'BorisMultipole':
        from . import dynamics_Boris_multipole as dynmul
        dynamics=dynmul.pusher_Boris_multipole(Dt=cc.Dt, N_sub_steps=cc.N_sub_steps, B_multip=cc.B_multip, B_skew=cc.B_skew)
    else:
        raise ValueError("""track_method should be 'Boris' or 'StrongBdip' or 'StrongBgen' or 'BorisMultipole'""")


    if cc.init_unif_flag==1:
        print("Adding inital %.2e electrons to the initial distribution" % cc.Nel_init_unif)
        MP_e.add_uniform_MP_distrib(cc.Nel_init_unif, cc.E_init_unif, cc.x_max_init_unif, cc.x_min_init_unif, cc.y_max_init_unif, cc.y_min_init_unif)


    if cc.init_unif_edens_flag==1:
        print("Adding inital %.2e electrons/m^3 to the initial distribution" % cc.init_unif_edens)
        MP_e.add_uniform_ele_density(n_ele=cc.init_unif_edens, E_init=cc.E_init_unif_edens,
                x_max=cc.x_max_init_unif_edens, x_min=cc.x_min_init_unif_edens,
                y_max=cc.y_max_init_unif_edens, y_min=cc.y_min_init_unif_edens)


    if cc.filename_init_MP_state!=-1 and cc.filename_init_MP_state is not None:
        print("Adding inital electrons from: %s" % cc.filename_init_MP_state)
        MP_e.add_from_file(cc.filename_init_MP_state)


    return (beamtim,
            MP_e,
            dynamics,impact_man,
            pyeclsaver,
            cc.gas_ion_flag,
            resgasion,
            cc.t_ion,
            spacech_ele,
            cc.t_sc_ON,
            cc.photoem_flag,
            phemiss,
            flag_presence_sec_beams,
            sec_beams_list,
            )


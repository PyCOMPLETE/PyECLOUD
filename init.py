#----------------------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                          PyECLOUD Version 6.3.0
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

from numpy import *
import beam_and_timing as beatim


from geom_impact_ellip import ellip_cham_geom_object

import sec_emission
from sec_emission_model_ECLOUD import SEY_model_ECLOUD
from sec_emission_model_accurate_low_ene import SEY_model_acc_low_ene
from sec_emission_model_ECLOUD_nunif import SEY_model_ECLOUD_non_unif
from sec_emission_model_cos_low_ener import SEY_model_cos_le
from sec_emission_model_flat_low_ener import SEY_model_flat_le
from sec_emission_model_from_file import SEY_model_from_file

import dynamics_dipole as dyndip
import dynamics_Boris_f2py as dynB
import dynamics_strong_B_generalized as dyngen

#import geom_impact_poly as gip


import MP_system as MPs
import space_charge_class as scc
import impact_management_class as imc
import pyecloud_saver as pysav
import gas_ionization_class as gic
import gen_photoemission_class as gpc

import parse_beam_file as pbf

qe=1.602176565e-19;
c=299792458.;

def read_parameter_files(pyecl_input_folder='./'):
    switch_model=0
    simulation_param_file='simulation_parameters.input'

    save_mp_state_time_file = -1

    stopfile = 'stop'

    dec_fact_out = 1

    init_unif_flag = 0
    Nel_init_unif = None
    E_init_unif = 0.
    x_max_init_unif = None
    x_min_init_unif = None
    y_max_init_unif = None
    y_min_init_unif = None

    chamb_type = 'ellip'
    filename_chm = None

    x_aper = None
    y_aper = None
    flag_detailed_MP_info=0
    flag_hist_impact_seg = 0

    track_method= 'StrongBdip'

    secondary_angle_distribution = 'cosine_2D'
    photoelectron_angle_distribution = 'cosine_2D'

    B = 0.   #Tesla (if B=-1 computed from energy and bending radius)
    bm_totlen= -1 #m


    B0x = 0.
    B0y = 0.
    B0z = 0.
    B_map_file = None
    Bz_map_file = None
    N_sub_steps = 1
    fact_Bmap = 1.
    B_zero_thrhld = None


    # photoemission parameters
    photoem_flag = 0
    inv_CDF_refl_photoem_file = -1
    k_pe_st = -1
    refl_frac = -1
    alimit= -1
    e_pe_sigma = -1
    e_pe_max = -1
    x0_refl = -1
    y0_refl = -1
    out_radius = -1

    # gas ionization parameters
    gas_ion_flag = 0
    P_nTorr=-1
    sigma_ion_MBarn=-1
    Temp_K=-1
    unif_frac=-1
    E_init_ion=-1

    N_mp_soft_regen = None
    N_mp_after_soft_regen = None
    Dx = 0.
    Dy = 0.
    betafx = None
    betafy = None


    flag_verbose_file=False
    flag_verbose_stdout=False

    secondary_beams_file_list = []

    phem_resc_fac = 0.9999

    dec_fac_secbeam_prof=1

    el_density_probes=[]

    save_simulation_state_time_file = -1

    # detailed histogram
    x_min_hist_det=None
    x_max_hist_det=None
    y_min_hist_det=None
    y_max_hist_det=None
    Dx_hist_det=None

    filename_init_MP_state = None

    sparse_solver = 'scipy_slu'

    B_multip = []
    B_skew = None

    PyPICmode = 'FiniteDifferences_ShortleyWeller'

    # uniform initial density
    init_unif_edens_flag = 0
    init_unif_edens = None
    E_init_unif_edens= None
    x_max_init_unif_edens = None
    x_min_init_unif_edens = None
    y_max_init_unif_edens = None
    y_min_init_unif_edens = None


    flag_assume_convex = True

    E0 = None
    s_param = None

    # multigrid parameters
    f_telescope = None
    target_grid = None
    N_nodes_discard = None
    N_min_Dh_main = None

    # secondary emission from file

    sey_file, range_extrapolate_right, delta_e, flag_factor_costheta, max_sey = [None]*5



    f=open(pyecl_input_folder+'/'+simulation_param_file)
    exec(f.read())
    f.close()


    f=open(pyecl_input_folder+'/'+machine_param_file)
    exec(f.read())
    f.close()

    f=open(pyecl_input_folder+'/'+secondary_emission_parameters_file)
    exec(f.read())
    f.close()

    b_par = pbf.beam_descr_from_fil(pyecl_input_folder+'/'+beam_parameters_file, betafx, Dx, betafy, Dy)

    flag_presence_sec_beams = False
    if len(secondary_beams_file_list)>0:
        flag_presence_sec_beams = True

    sec_b_par_list=[]
    if flag_presence_sec_beams:
        for sec_b_file in secondary_beams_file_list:
            sec_b_par_list.append(pbf.beam_descr_from_fil(pyecl_input_folder+'/'+sec_b_file, betafx, Dx, betafy, Dy))

    if B==-1:
        B = 2*pi*b_par.beta_rel*b_par.energy_J/(c*qe*bm_totlen)

    filen_main_outp = 'Pyecltest'

    return (
        b_par,
        x_aper,
        y_aper,
        B,
        gas_ion_flag,
        P_nTorr,
        sigma_ion_MBarn,
        Temp_K,
        unif_frac,
        E_init_ion,
        Emax,
        del_max,
        R0, E_th,
        sigmafit,
        mufit,
        Dt,
        t_end,
        lam_th,
        t_ion,
        N_mp_max,
        N_mp_regen,
        N_mp_after_regen,
        fact_split,
        fact_clean,
        nel_mp_ref_0,
        Nx_regen,
        Ny_regen,
        Nvx_regen,
        Nvy_regen,
        Nvz_regen,
        regen_hist_cut,
        N_mp_regen_low,
        Dt_sc,
        Dh_sc,
        t_sc_ON,
        Dx_hist,
        r_center,
        scrub_en_th,
        progress_path,
        logfile_path,
        flag_movie,
        flag_sc_movie,
        Dt_En_hist,
        Nbin_En_hist,
        En_hist_max,
        photoem_flag,
        inv_CDF_refl_photoem_file,
        k_pe_st,
        refl_frac,
        alimit,
        e_pe_sigma,
        e_pe_max,
        x0_refl,
        y0_refl,
        out_radius,
        switch_model,
        switch_no_increase_energy,
        thresh_low_energy,
        save_mp_state_time_file,
        init_unif_flag,
        Nel_init_unif,
        E_init_unif,
        x_max_init_unif,
        x_min_init_unif,
        y_max_init_unif,
        y_min_init_unif,
        chamb_type,
        filename_chm,
        flag_detailed_MP_info,
        flag_hist_impact_seg,
        track_method,
        secondary_angle_distribution,
        photoelectron_angle_distribution,
        B0x,
        B0y,
        B0z,
        B_map_file,
        Bz_map_file,
        N_sub_steps,
        fact_Bmap,
        B_zero_thrhld,
        N_mp_soft_regen,
        N_mp_after_soft_regen,
        flag_verbose_file,
        flag_verbose_stdout,
        flag_presence_sec_beams,
        sec_b_par_list,
        phem_resc_fac,
        dec_fac_secbeam_prof,
        el_density_probes,
        save_simulation_state_time_file,
        x_min_hist_det,
        x_max_hist_det,
        y_min_hist_det,
        y_max_hist_det,
        Dx_hist_det,
        dec_fact_out,
        stopfile,
        sparse_solver,
        B_multip,
        B_skew,
        PyPICmode,
        filename_init_MP_state,
        init_unif_edens_flag,
        init_unif_edens,
        E_init_unif_edens,
        x_max_init_unif_edens,
        x_min_init_unif_edens,
        y_max_init_unif_edens,
        y_min_init_unif_edens,
        flag_assume_convex,
        E0,
        s_param,
        filen_main_outp,
        f_telescope,
        target_grid,
        N_nodes_discard,
        N_min_Dh_main,
        sey_file,
        range_extrapolate_right,
        delta_e,
        flag_factor_costheta,
        max_sey,
    )


def read_input_files_and_init_components(pyecl_input_folder='./', **kwargs):
    (
        b_par,
        x_aper,
        y_aper,
        B,
        gas_ion_flag,
        P_nTorr,
        sigma_ion_MBarn,
        Temp_K,
        unif_frac,
        E_init_ion,
        Emax,
        del_max,
        R0, E_th,
        sigmafit,
        mufit,
        Dt,
        t_end,
        lam_th,
        t_ion,
        N_mp_max,
        N_mp_regen,
        N_mp_after_regen,
        fact_split,
        fact_clean,
        nel_mp_ref_0,
        Nx_regen,
        Ny_regen,
        Nvx_regen,
        Nvy_regen,
        Nvz_regen,
        regen_hist_cut,
        N_mp_regen_low,
        Dt_sc,
        Dh_sc,
        t_sc_ON,
        Dx_hist,
        r_center,
        scrub_en_th,
        progress_path,
        logfile_path,
        flag_movie,
        flag_sc_movie,
        Dt_En_hist,
        Nbin_En_hist,
        En_hist_max,
        photoem_flag,
        inv_CDF_refl_photoem_file,
        k_pe_st,
        refl_frac,
        alimit,
        e_pe_sigma,
        e_pe_max,
        x0_refl,
        y0_refl,
        out_radius,
        switch_model,
        switch_no_increase_energy,
        thresh_low_energy,
        save_mp_state_time_file,
        init_unif_flag,
        Nel_init_unif,
        E_init_unif,
        x_max_init_unif,
        x_min_init_unif,
        y_max_init_unif,
        y_min_init_unif,
        chamb_type,
        filename_chm,
        flag_detailed_MP_info,
        flag_hist_impact_seg,
        track_method,
        secondary_angle_distribution,
        photoelectron_angle_distribution,
        B0x,
        B0y,
        B0z,
        B_map_file,
        Bz_map_file,
        N_sub_steps,
        fact_Bmap,
        B_zero_thrhld,
        N_mp_soft_regen,
        N_mp_after_soft_regen,
        flag_verbose_file,
        flag_verbose_stdout,
        flag_presence_sec_beams,
        sec_b_par_list,
        phem_resc_fac,
        dec_fac_secbeam_prof,
        el_density_probes,
        save_simulation_state_time_file,
        x_min_hist_det,
        x_max_hist_det,
        y_min_hist_det,
        y_max_hist_det,
        Dx_hist_det,
        dec_fact_out,
        stopfile,
        sparse_solver,
        B_multip,
        B_skew,
        PyPICmode,
        filename_init_MP_state,
        init_unif_edens_flag,
        init_unif_edens,
        E_init_unif_edens,
        x_max_init_unif_edens,
        x_min_init_unif_edens,
        y_max_init_unif_edens,
        y_min_init_unif_edens,
        flag_assume_convex,
        E0,
        s_param,
        filen_main_outp,
        f_telescope,
        target_grid,
        N_nodes_discard,
        N_min_Dh_main,
        sey_file,
        range_extrapolate_right,
        delta_e,
        flag_factor_costheta,
        max_sey,
        ) = read_parameter_files(pyecl_input_folder)



    for attr in kwargs.keys():
            print 'Ecloud init. From kwargs: %s = %s'%(attr, repr(kwargs[attr]))
            tmpattr = kwargs[attr]
            exec('%s=tmpattr'%attr)





    ##########################################

    pyeclsaver=pysav.pyecloud_saver(logfile_path)

    if switch_model=='ECLOUD_nunif':
        flag_non_unif_sey = 1
    else:
        flag_non_unif_sey = 0

    if chamb_type=='ellip':
        chamb=ellip_cham_geom_object(x_aper, y_aper, flag_verbose_file=flag_verbose_file)
    elif chamb_type=='polyg' or chamb_type=='polyg_cython':

        import geom_impact_poly_fast_impact as gipfi
        chamb=gipfi.polyg_cham_geom_object(filename_chm, flag_non_unif_sey,
                                     flag_verbose_file=flag_verbose_file, flag_verbose_stdout=flag_verbose_stdout, flag_assume_convex=flag_assume_convex)
    elif chamb_type=='polyg_numpy':
        raise ValueError("chamb_type='polyg_numpy' not supported anymore")
        #~ chamb=gip.polyg_cham_geom_object(filename_chm, flag_non_unif_sey,
                         #~ flag_verbose_file=flag_verbose_file, flag_verbose_stdout=flag_verbose_stdout)
    elif chamb_type=='rect':
        import geom_impact_rect_fast_impact as girfi
        chamb =  girfi.rect_cham_geom_object(x_aper, y_aper, flag_verbose_file=flag_verbose_file, flag_verbose_stdout=flag_verbose_stdout)
    else:
        raise ValueError('Chamber type not recognized (choose: ellip/rect/polyg)')


    MP_e=MPs.MP_system(N_mp_max, nel_mp_ref_0, fact_split, fact_clean,
                       N_mp_regen_low, N_mp_regen, N_mp_after_regen,
                       Dx_hist, Nx_regen, Ny_regen, Nvx_regen, Nvy_regen, Nvz_regen, regen_hist_cut, chamb,
                       N_mp_soft_regen=N_mp_soft_regen, N_mp_after_soft_regen=N_mp_after_soft_regen)

    beamtim=beatim.beam_and_timing(b_par.flag_bunched_beam, b_par.fact_beam, b_par.coast_dens, b_par.beam_field_file,lam_th,
                 b_spac=b_par.b_spac, sigmaz=b_par.sigmaz,t_offs=b_par.t_offs, filling_pattern_file=b_par.filling_pattern_file, Dt=Dt, t_end=t_end,
                 beam_long_prof_file=b_par.beam_long_prof_file, Dh_beam_field=b_par.Dh_beam_field, f_telescope_beam = b_par.f_telescope_beam,
                 target_grid_beam = b_par.target_grid_beam, N_nodes_discard_beam = b_par.N_nodes_discard_beam, N_min_Dh_main_beam = b_par.N_min_Dh_main_beam,
                 chamb=chamb,  sigmax=b_par.sigmax, sigmay=b_par.sigmay,
                 x_beam_pos = b_par.x_beam_pos, y_beam_pos = b_par.y_beam_pos, save_beam_field_file_as=b_par.save_beam_field_file_as,
                 Nx=b_par.Nx, Ny=b_par.Ny, nimag=b_par.nimag, progress_mapgen_file = (progress_path+'_mapgen'))

    if sparse_solver=='klu':
        print '''sparse_solver: 'klu' no longer supported --> going to PyKLU'''
        sparse_solver='PyKLU'

    spacech_ele = scc.space_charge(chamb, Dh_sc, Dt_sc=Dt_sc, sparse_solver=sparse_solver, PyPICmode=PyPICmode,
                        f_telescope = f_telescope, target_grid = target_grid, N_nodes_discard = N_nodes_discard, N_min_Dh_main = N_min_Dh_main)

    sec_beams_list=[]
    if flag_presence_sec_beams:
        N_sec_beams = len(sec_b_par_list)
        for ii in xrange(N_sec_beams):
            print 'Initialize secondary beam %d/%d'%(ii+1, N_sec_beams)
            sb_par = sec_b_par_list[ii]
            sec_beams_list.append(beatim.beam_and_timing(sb_par.flag_bunched_beam, sb_par.fact_beam, sb_par.coast_dens, sb_par.beam_field_file,lam_th,
                 b_spac=sb_par.b_spac, sigmaz=sb_par.sigmaz,t_offs=sb_par.t_offs, filling_pattern_file=sb_par.filling_pattern_file, Dt=Dt, t_end=t_end,
                 beam_long_prof_file=sb_par.beam_long_prof_file, Dh_beam_field=sb_par.Dh_beam_field, f_telescope_beam = sb_par.f_telescope_beam,
                 target_grid_beam = sb_par.target_grid_beam, N_nodes_discard_beam = sb_par.N_nodes_discard_beam, N_min_Dh_main_beam = sb_par.N_min_Dh_main_beam,
                 chamb=chamb,  sigmax=sb_par.sigmax, sigmay=sb_par.sigmay,
                 x_beam_pos = sb_par.x_beam_pos, y_beam_pos = sb_par.y_beam_pos, save_beam_field_file_as=sb_par.save_beam_field_file_as,
                 flag_secodary_beam = True, t_primary_beam = beamtim.t,
                 Nx=sb_par.Nx, Ny=sb_par.Ny, nimag=sb_par.nimag, progress_mapgen_file = (progress_path+('_mapgen_sec_%d'%ii))))

    kwargs = {}

    if E0 is not None:
        kwargs.update({'E0':E0})
        #If E0 is not provided use default value for each object

    if s_param is not None:
        if switch_model==0 or switch_model=='ECLOUD':
            kwargs.update({'s':s_param})
        else:
            raise ValueError('s parameter can be changed only in the ECLOUD sec. emission model!')

    if switch_model==0 or switch_model=='ECLOUD':
        sey_mod=SEY_model_ECLOUD(Emax,del_max,R0,**kwargs)
    elif switch_model==1 or switch_model=='ACC_LOW':
        sey_mod=SEY_model_acc_low_ene(Emax,del_max,R0,**kwargs)
    elif switch_model=='ECLOUD_nunif':
        sey_mod=SEY_model_ECLOUD_non_unif(chamb, Emax,del_max,R0,**kwargs)
    elif switch_model=='cos_low_ene':
        sey_mod=SEY_model_cos_le(Emax,del_max,R0,**kwargs)
    elif switch_model=='flat_low_ene':
        sey_mod=SEY_model_flat_le(Emax,del_max,R0)
    elif switch_model == 'from_file':
        sey_mod = SEY_model_from_file(sey_file, range_extrapolate_right, delta_e, flag_factor_costheta, max_sey)


    secondary_angle_dist_func = {
        'cosine_3D': sec_emission.velocities_angle_cosine_3D,
        'cosine_2D': sec_emission.velocities_angle_cosine_2D,
    }[secondary_angle_distribution]

    photoelectron_angle_dist_func = {
        'cosine_3D': sec_emission.velocities_angle_cosine_3D,
        'cosine_2D': sec_emission.velocities_angle_cosine_2D,
    }[photoelectron_angle_distribution]



    flag_seg = (flag_hist_impact_seg==1)

    impact_man=imc.impact_management(switch_no_increase_energy, chamb, sey_mod, E_th, sigmafit, mufit,
                 Dx_hist, scrub_en_th, Nbin_En_hist, En_hist_max, thresh_low_energy=thresh_low_energy,
                 flag_seg=flag_seg, angle_dist_func=secondary_angle_dist_func)


    #resgasion_sec_beam_list=[]
    if gas_ion_flag==1:
        resgasion=gic.residual_gas_ionization(unif_frac, P_nTorr, sigma_ion_MBarn,Temp_K,chamb,E_init_ion)
    else:
        resgasion=None



    if photoem_flag==1:
        phemiss=gpc.photoemission(inv_CDF_refl_photoem_file, k_pe_st, refl_frac, e_pe_sigma, e_pe_max,alimit, \
                x0_refl, y0_refl, out_radius, chamb, phem_resc_fac, photoelectron_angle_dist_func)
    else:
        phemiss=None

    pyeclsaver.start_observing(MP_e, beamtim, impact_man,
                 r_center, Dt_En_hist, logfile_path, progress_path, flag_detailed_MP_info=flag_detailed_MP_info,
                 flag_movie=flag_movie, flag_sc_movie=flag_sc_movie, save_mp_state_time_file=save_mp_state_time_file,
                 flag_presence_sec_beams=flag_presence_sec_beams, sec_beams_list=sec_beams_list, dec_fac_secbeam_prof=dec_fac_secbeam_prof,
                 el_density_probes=el_density_probes, save_simulation_state_time_file = save_simulation_state_time_file,
                 x_min_hist_det=x_min_hist_det, x_max_hist_det=x_max_hist_det, y_min_hist_det=y_min_hist_det, y_max_hist_det=y_max_hist_det,
                 Dx_hist_det=Dx_hist_det, dec_fact_out=dec_fact_out, stopfile=stopfile, filen_main_outp=filen_main_outp)






    if track_method == 'Boris':
        dynamics=dynB.pusher_Boris(Dt, B0x, B0y, B0z, \
                 B_map_file, fact_Bmap,  Bz_map_file,N_sub_steps=N_sub_steps)
    elif track_method == 'StrongBdip':
        dynamics=dyndip.pusher_dipole_magnet(Dt,B)
    elif track_method == 'StrongBgen':
        dynamics=dyngen.pusher_strong_B_generalized(Dt, B0x, B0y,  \
                 B_map_file, fact_Bmap, B_zero_thrhld)
    elif track_method == 'BorisMultipole':
        import dynamics_Boris_multipole as dynmul
        dynamics=dynmul.pusher_Boris_multipole(Dt=Dt, N_sub_steps=N_sub_steps, B_multip=B_multip, B_skew=B_skew)
    else:
        raise ValueError("""track_method should be 'Boris' or 'StrongBdip' or 'StrongBgen' or 'BorisMultipole'""")



    if init_unif_flag==1:
        print "Adding inital %.2e electrons to the initial distribution"%Nel_init_unif
        MP_e.add_uniform_MP_distrib(Nel_init_unif, E_init_unif, x_max_init_unif, x_min_init_unif, y_max_init_unif, y_min_init_unif)


    if init_unif_edens_flag==1:
        print "Adding inital %.2e electrons/m^3 to the initial distribution"%init_unif_edens
        MP_e.add_uniform_ele_density(n_ele=init_unif_edens, E_init=E_init_unif_edens,
        x_max=x_max_init_unif_edens, x_min=x_min_init_unif_edens,
        y_max=y_max_init_unif_edens, y_min=y_min_init_unif_edens)


    if filename_init_MP_state!=-1 and filename_init_MP_state is not None:
        print "Adding inital electrons from: %s"%filename_init_MP_state
        MP_e.add_from_file(filename_init_MP_state)


    return beamtim,MP_e, dynamics,impact_man, pyeclsaver, \
        gas_ion_flag, resgasion, t_ion, \
        spacech_ele,t_sc_ON, photoem_flag, phemiss,\
        flag_presence_sec_beams, sec_beams_list

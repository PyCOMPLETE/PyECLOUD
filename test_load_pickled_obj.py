import pickle
from . import pyecloud_saver as pysav

filename_sim_state = 'simulation_state_0.pkl'

r_center = 1e-3
Dt_En_hist = 25e-9
logfile_path= 'logfile_restarted.txt'
progress_path= 'progress_restarted.txt'
flag_detailed_MP_info = 1
flag_movie = 0
flag_sc_movie=  0
save_mp_state_time_file = -1
dec_fac_secbeam_prof=1





with open(filename_sim_state, 'rb') as fid:
    dict_state = pickle.load(fid)


beamtim=dict_state['beamtim']
MP_e=dict_state['MP_e']
dynamics=dict_state['dynamics']
impact_man=dict_state['impact_man']
gas_ion_flag=dict_state['gas_ion_flag']
resgasion=dict_state['resgasion']
t_ion=dict_state['t_ion']
spacech_ele=dict_state['spacech_ele']
t_sc_ON=dict_state['t_sc_ON']
photoem_flag=dict_state['photoem_flag']
phemiss=dict_state['phemiss']
flag_presence_sec_beams=dict_state['flag_presence_sec_beams']
sec_beams_list=dict_state['sec_beams_list']
el_density_probes = []







pyeclsaver=pysav.pyecloud_saver(MP_e, beamtim, impact_man,
                 r_center, Dt_En_hist, logfile_path, progress_path, flag_detailed_MP_info=flag_detailed_MP_info,
                 flag_movie=flag_movie, flag_sc_movie=flag_sc_movie, save_mp_state_time_file=save_mp_state_time_file,
                 flag_presence_sec_beams=flag_presence_sec_beams, sec_beams_list=sec_beams_list, dec_fac_secbeam_prof=dec_fac_secbeam_prof,
                 el_density_probes=el_density_probes)

print('Start timestep iter')

## simulation
while not beamtim.end_simulation():

    beamtim.next_time_step()

    if flag_presence_sec_beams:
        for sec_beam in sec_beams_list:
            sec_beam.next_time_step()


    ## compute beam electric field (main and secondary beams)
    Ex_n_beam, Ey_n_beam = beamtim.get_beam_eletric_field(MP_e)

    if flag_presence_sec_beams:
        for sec_beam in sec_beams_list:
            Ex_n_secbeam, Ey_n_secbeam = sec_beam.get_beam_eletric_field(MP_e)
            Ex_n_beam+=Ex_n_secbeam
            Ey_n_beam+=Ey_n_secbeam

    ## compute electron space charge electric field
    Ex_sc_n, Ey_sc_n = spacech_ele.get_sc_eletric_field(MP_e)

    ## Total electric field
    Ex_n=Ex_sc_n+Ex_n_beam;
    Ey_n=Ey_sc_n+Ey_n_beam;

    ## save position before motion step
    old_pos=MP_e.get_positions()

    ## motion
    MP_e = dynamics.step(MP_e, Ex_n, Ey_n);

    ## impacts: backtracking and secondary emission
    MP_e = impact_man.backtrack_and_second_emiss(old_pos, MP_e)


    ## gas ionization (main and secondary beams)
    if(beamtim.tt_curr<t_ion and gas_ion_flag==1):
        MP_e = resgasion.generate(MP_e, beamtim.lam_t_curr, beamtim.Dt,beamtim.sigmax, beamtim.sigmay,
                                x_beam_pos = beamtim.x_beam_pos, y_beam_pos = beamtim.y_beam_pos)
        if flag_presence_sec_beams:
            for sec_beam in sec_beams_list:
                MP_e = resgasion.generate(MP_e, sec_beam.lam_t_curr, sec_beam.Dt,sec_beam.sigmax, sec_beam.sigmay,
                                x_beam_pos = sec_beam.x_beam_pos, y_beam_pos = sec_beam.y_beam_pos)



    ## photoemission (main and secondary beams)
    if(photoem_flag==1):
        lam_curr_phem = beamtim.lam_t_curr
        if flag_presence_sec_beams:
            for sec_beam in sec_beams_list:
                lam_curr_phem += sec_beam.lam_t_curr
        MP_e = phemiss.generate(MP_e, lam_curr_phem, beamtim.Dt)


    ## Compute space charge field
    if (beamtim.tt_curr>t_sc_ON):
        spacech_ele.recompute_spchg_efield(MP_e, t_curr=beamtim.tt_curr)


    ## savings
    impact_man = pyeclsaver.witness(MP_e, beamtim, spacech_ele, impact_man, dynamics,gas_ion_flag,
                                    resgasion,t_ion,t_sc_ON, photoem_flag, phemiss,
                                    flag_presence_sec_beams,sec_beams_list)


    ## every bunch passage
    if beamtim.flag_new_bunch_pass:

        ## clean
        MP_e.clean_small_MPs()

        ## regeneration
        MP_e.check_for_regeneration()

        ## soft regeneration
        MP_e.check_for_soft_regeneration()

        print('**** Done pass_numb = %d/%d\n'%(beamtim.pass_numb,beamtim.N_pass_tot))




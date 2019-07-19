#!/afs/cern.ch/project/uslarp/opt/lxplus64/Python-2.7.2/bin/python

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


import init as init
import cPickle
import numpy as np
import os


class BuildupSimulation(object):

    def __init__(self, pyecl_input_folder='./', skip_beam=False, skip_spacech_ele=False,
                 skip_pyeclsaver=False, ignore_kwargs=[], spacech_ele=None, **kwargs):

        print 'PyECLOUD Version 8.1.0'
        beamtim, spacech_ele, t_sc_ON, flag_presence_sec_beams, sec_beams_list, \
            config_dict, flag_multiple_clouds, cloud_list, checkpoint_folder, \
            cross_ion = init.read_input_files_and_init_components(\
                pyecl_input_folder=pyecl_input_folder,
                skip_beam=skip_beam,
                skip_pyeclsaver=skip_pyeclsaver,
                skip_spacech_ele=skip_spacech_ele,
                spacech_ele=spacech_ele,
                ignore_kwargs=ignore_kwargs,
                **kwargs)

        self.config_dict = config_dict
        self.beamtim = beamtim
        self.spacech_ele = spacech_ele
        self.t_sc_ON = t_sc_ON
        self.flag_presence_sec_beams = flag_presence_sec_beams
        self.sec_beams_list = sec_beams_list
        self.flag_multiple_clouds = flag_multiple_clouds
        self.cloud_list = cloud_list
        self.chamb = cloud_list[0].impact_man.chamb
        self.checkpoint_folder = checkpoint_folder
        self.cross_ion = cross_ion

        # Checking if there are saved checkpoints
        if self.checkpoint_folder is not None:
            if os.path.isdir(self.checkpoint_folder):
                if len(os.listdir(self.checkpoint_folder)) == 1:
                    print('Loading from checkpoint...')
                    # Selecting latest checkpoint
                    checkpoint = os.listdir(self.checkpoint_folder)[0]
                    self.load_checkpoint(filename_simulation_checkpoint=checkpoint, load_from_folder=self.checkpoint_folder)
                elif len(os.listdir(self.checkpoint_folder)) == 0:
                    print('No checkpoint found, starting new simulation...')
                else:
                    raise ValueError('More than one checkpoint found in %s'%self.checkpoint_folder)


    def run(self, t_end_sim=None):

        beamtim = self.beamtim

        flag_presence_sec_beams = self.flag_presence_sec_beams
        sec_beams_list = self.sec_beams_list

        print 'Start timestep iter'

        ## simulation
        while not beamtim.end_simulation():

            if t_end_sim is not None and beamtim.tt_curr is not None:
                if beamtim.tt_curr >= t_end_sim:
                    print 'Reached user defined t_end_sim --> Ending simulation'
                    break

            beamtim.next_time_step()

            if flag_presence_sec_beams:
                for sec_beam in sec_beams_list:
                    sec_beam.next_time_step()

            self.sim_time_step()

            if beamtim.flag_new_bunch_pass:
                print '**** Done pass_numb = %d/%d\n'%(beamtim.pass_numb, beamtim.N_pass_tot)


    def sim_time_step(self, beamtim_obj=None, Dt_substep_custom=None, N_sub_steps_custom=None, kick_mode_for_beam_field=False,
                      force_recompute_space_charge=False, skip_MP_cleaning=False, skip_MP_regen=False):

        if beamtim_obj is not None:
            beamtim = beamtim_obj
        else:
            beamtim = self.beamtim

        spacech_ele = self.spacech_ele
        t_sc_ON = self.t_sc_ON
        flag_presence_sec_beams = self.flag_presence_sec_beams
        sec_beams_list = self.sec_beams_list
        flag_multiple_clouds = self.flag_multiple_clouds
        cloud_list = self.cloud_list
        cross_ion = self.cross_ion

        flag_recompute_space_charge = spacech_ele.check_for_recomputation(t_curr=beamtim.tt_curr)

        # Loop over clouds: gather fields, move, generate new MPs
        for i_cloud, cloud in enumerate(cloud_list):
            MP_e = cloud.MP_e
            dynamics = cloud.dynamics
            impact_man = cloud.impact_man
            gas_ion_flag = cloud.gas_ion_flag
            resgasion = cloud.resgasion
            t_ion = cloud.t_ion
            photoem_flag = cloud.photoem_flag
            phemiss = cloud.phemiss

            ## Compute beam electric field (main and secondary beams)
            Ex_n_beam, Ey_n_beam = beamtim.get_beam_eletric_field(MP_e)

            if flag_presence_sec_beams:
                for sec_beam in sec_beams_list:
                    Ex_n_secbeam, Ey_n_secbeam = sec_beam.get_beam_eletric_field(MP_e)
                    Ex_n_beam += Ex_n_secbeam
                    Ey_n_beam += Ey_n_secbeam

            ## Compute electron space charge electric field
            Ex_sc_n, Ey_sc_n = spacech_ele.get_sc_eletric_field(MP_e)

            if kick_mode_for_beam_field:
                if Dt_substep_custom is None or N_sub_steps_custom is None:
                    raise ValueError("""Kick mode can be used only with custom time steps!""")
                MP_e.vx_mp[:MP_e.N_mp] += Ex_n_beam * (Dt_substep_custom * N_sub_steps_custom) * MP_e.charge / MP_e.mass
                MP_e.vy_mp[:MP_e.N_mp] += Ey_n_beam * (Dt_substep_custom * N_sub_steps_custom) * MP_e.charge / MP_e.mass
                # Electric field for dynamics step
                Ex_n = Ex_sc_n
                Ey_n = Ey_sc_n
            else:
                # Electric field for dynamics step
                Ex_n = Ex_sc_n + Ex_n_beam
                Ey_n = Ey_sc_n + Ey_n_beam

            ## Save position before motion step
            old_pos = MP_e.get_positions()

            ## Motion
            if Dt_substep_custom is None and N_sub_steps_custom is None and beamtim.flag_unif_Dt:
                # Standard simulation mode
                MP_e = dynamics.step(MP_e, Ex_n, Ey_n)
            elif Dt_substep_custom is None and N_sub_steps_custom is None and not(beamtim.flag_unif_Dt):
                # Dt from non-uniform beam profile
                if self.config_dict['track_method'] not in ['Boris', 'BorisMultipole']:
                    raise ValueError("""track_method should be 'Boris' or 'BorisMultipole' to use custom substeps!""")
                Dt_substep_target = cloud.dynamics.Dt / cloud.dynamics.N_sub_steps
                N_substeps_curr = np.round(beamtim.Dt_curr / Dt_substep_target)
                Dt_substep_curr = beamtim.Dt_curr / N_substeps_curr
                MP_e = dynamics.stepcustomDt(MP_e, Ex_n, Ey_n, Dt_substep=Dt_substep_curr, N_sub_steps=N_substeps_curr)
            else:
                # Custom steps and substeps provided as arguments of sim_time_step
                if self.config_dict['track_method'] not in ['Boris', 'BorisMultipole']:
                    raise ValueError("""track_method should be 'Boris' or 'BorisMultipole' to use custom substeps!""")
                MP_e = dynamics.stepcustomDt(MP_e, Ex_n, Ey_n, Dt_substep=Dt_substep_custom, N_sub_steps=N_sub_steps_custom)

            ## Impacts: backtracking and secondary emission
            MP_e = impact_man.backtrack_and_second_emiss(old_pos, MP_e, beamtim.tt_curr)

            ## Evolve SEY module (e.g. charge decay for insulators
            impact_man.sey_mod.SEY_model_evol(Dt=beamtim.Dt_curr)

            ## Gas ionization (main and secondary beams)
            if(beamtim.tt_curr < t_ion and gas_ion_flag == 1):
                MP_e = resgasion.generate(MP_e, beamtim.lam_t_curr, beamtim.Dt_curr, beamtim.sigmax, beamtim.sigmay,
                                          x_beam_pos=beamtim.x_beam_pos, y_beam_pos=beamtim.y_beam_pos)
                if flag_presence_sec_beams:
                    for sec_beam in sec_beams_list:
                        MP_e = resgasion.generate(MP_e, sec_beam.lam_t_curr, sec_beam.Dt_curr, sec_beam.sigmax, sec_beam.sigmay,
                                                  x_beam_pos=sec_beam.x_beam_pos, y_beam_pos=sec_beam.y_beam_pos)

            ## Photoemission (main and secondary beams)
            if (photoem_flag != 0):
                lam_curr_phem = beamtim.lam_t_curr
                if flag_presence_sec_beams:
                    for sec_beam in sec_beams_list:
                        lam_curr_phem += sec_beam.lam_t_curr
                phemiss.generate(MP_e, lam_curr_phem, beamtim.Dt_curr)


        ## Cross_ionization
        if cross_ion is not None:
           cross_ion.generate(Dt=beamtim.Dt_curr, cloud_list=cloud_list)


        ## Compute space charge field
        for cloud in cloud_list:
            if ((beamtim.tt_curr > t_sc_ON) and flag_recompute_space_charge) or force_recompute_space_charge:
                flag_reset = cloud is cloud_list[0] # The first cloud resets the distribution
                flag_solve = cloud is cloud_list[-1] # The last cloud computes the fields
                spacech_ele.recompute_spchg_efield(MP_e, flag_solve=flag_solve, flag_reset=flag_reset)

                # Copy rho to cloud
                cloud.rho = spacech_ele.rho - sum([cl.rho for cl in cloud_list[:i_cloud]])

        ## Saving output
        # We want to save and clean MP only after iteration on all clouds is completed
        # (e.g. to have consistent space charge state)
        for cloud in cloud_list:
            if cloud.pyeclsaver is not None:
                # if Dt_substep_custom is not None or N_sub_steps_custom is not None:
                #     raise ValueError('Saving with custom steps not implemented!')
                cloud.impact_man = cloud.pyeclsaver.witness(cloud.MP_e, beamtim, spacech_ele, cloud.impact_man, cloud.dynamics,
                                                            cloud.gas_ion_flag, cloud.resgasion, cloud.t_ion, t_sc_ON,
                                                            cloud.photoem_flag, cloud.phemiss, flag_presence_sec_beams,
                                                            sec_beams_list, cloud_list, buildup_sim=self, cross_ion=cross_ion,
                                                            rho_cloud=cloud.rho)

        ## Cleaning and regeneration
        for cloud in cloud_list:
            ## Every bunch passage
            if beamtim.flag_new_bunch_pass:
                ## Clean
                if not skip_MP_cleaning:
                    cloud.MP_e.clean_small_MPs(cloud.name)
                if not skip_MP_regen:
                    ## Regeneration
                    cloud.MP_e.check_for_regeneration(cloud.name)
                    ## Soft regeneration
                    cloud.MP_e.check_for_soft_regeneration(cloud.name)

            cloud.MP_e.check_for_async_regeneration(cloud.name)


    def load_state(self, filename_simulation_state, force_disable_save_simulation_state=True, filen_main_outp='Pyecltest_restarted', load_from_folder='./'):  # , reset_pyeclsaver = True):

        with open(load_from_folder + filename_simulation_state, 'rb') as fid:
            dict_state = cPickle.load(fid)

        self.beamtim = dict_state['beamtim']
        self.spacech_ele = dict_state['spacech_ele']
        self.t_sc_ON = dict_state['t_sc_ON']
        self.flag_presence_sec_beams = dict_state['flag_presence_sec_beams']
        self.sec_beams_list = dict_state['sec_beams_list']

        self.flag_multiple_clouds = dict_state['flag_multiple_clouds']

        for i_cloud, new_cloud in enumerate(self.cloud_list):
            new_pyeclsaver = new_cloud.pyeclsaver
            self.cloud_list[i_cloud] = dict_state['cloud_list'][i_cloud] # Replace new_cloud with saved cloud
            cloud = self.cloud_list[i_cloud]

            # if reset_pyeclsaver or cloud.pyeclsaver is None:
            cloud.pyeclsaver = new_pyeclsaver

            if force_disable_save_simulation_state:
                cloud.pyeclsaver.flag_save_simulation_state = False

            if filen_main_outp is not None:
                filen_outp_ext = filen_main_outp.split('Pyecltest')[-1]
                filen_outp_root = cloud.pyeclsaver.filen_main_outp.split('.mat')[0]
                cloud.pyeclsaver.filen_main_outp = filen_outp_root + filen_outp_ext + '.mat'

        print 'Restoring PyPIC LU object...'
        self.spacech_ele.PyPICobj.build_sparse_solver()
        print 'Done reload.'
        return dict_state


    def load_checkpoint(self, filename_simulation_checkpoint, load_from_folder='./'):
        print('Reloading from checkpoint: %s...' % (load_from_folder + filename_simulation_checkpoint))

        i_checkp = int(filename_simulation_checkpoint.split('.pkl')[0].split('_')[-1])
        dict_state = self.load_state(filename_simulation_checkpoint, force_disable_save_simulation_state=False, filen_main_outp=None, load_from_folder=load_from_folder)

        for cloud in self.cloud_list:
            cloud.pyeclsaver.load_from_output(last_t=self.beamtim.tt_curr)
            cloud.pyeclsaver.i_checkp = i_checkp + 1
            cloud.pyeclsaver.t_last_checkp = self.beamtim.tt_curr
            cloud.pyeclsaver.t_last_En_hist = dict_state['t_last_En_hist']

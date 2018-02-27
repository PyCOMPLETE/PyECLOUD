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
#                   PyECLOUD Version 6.7.2
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


class BuildupSimulation(object):
    def __init__(self, pyecl_input_folder='./', **kwargs):


        print 'PyECLOUD Version 6.7.2'
        beamtim, spacech_ele, t_sc_ON, flag_presence_sec_beams, sec_beams_list, \
        config_dict, flag_multiple_clouds, cloud_list = init.read_input_files_and_init_components(\
                                                    pyecl_input_folder=pyecl_input_folder, **kwargs)

        
        self.beamtim = beamtim
        self.spacech_ele = spacech_ele
        self.t_sc_ON = t_sc_ON
        self.flag_presence_sec_beams = flag_presence_sec_beams
        self.sec_beams_list = sec_beams_list
        self.flag_multiple_clouds = flag_multiple_clouds
        self.cloud_list = cloud_list

    def run(self, t_end_sim = None):

        beamtim = self.beamtim
        spacech_ele = self.spacech_ele
        t_sc_ON = self.t_sc_ON
        flag_presence_sec_beams = self.flag_presence_sec_beams
        sec_beams_list = self.sec_beams_list
        flag_multiple_clouds = self.flag_multiple_clouds
        cloud_list = self.cloud_list

        print 'Start timestep iter'

        ## simulation
        while not beamtim.end_simulation():

            beamtim.next_time_step()

            if flag_presence_sec_beams:
                for sec_beam in sec_beams_list:
                    sec_beam.next_time_step()


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
                    MP_e = resgasion.generate(MP_e, beamtim.lam_t_curr, beamtim.Dt, beamtim.sigmax, beamtim.sigmay,
                                              x_beam_pos = beamtim.x_beam_pos, y_beam_pos = beamtim.y_beam_pos)
                    if flag_presence_sec_beams:
                        for sec_beam in sec_beams_list:
                            MP_e = resgasion.generate(MP_e, sec_beam.lam_t_curr, sec_beam.Dt, sec_beam.sigmax, sec_beam.sigmay,
                                                      x_beam_pos = sec_beam.x_beam_pos, y_beam_pos = sec_beam.y_beam_pos)

                ## photoemission (main and secondary beams)
                if (photoem_flag != 0):
                    lam_curr_phem = beamtim.lam_t_curr
                    if flag_presence_sec_beams:
                        for sec_beam in sec_beams_list:
                            lam_curr_phem += sec_beam.lam_t_curr
                    phemiss.generate(MP_e, lam_curr_phem, beamtim.Dt)

                # Compute space charge field
                if (beamtim.tt_curr>t_sc_ON):
                    flag_reset = cloud is cloud_list[0] # The first cloud resets the distribution
                    flag_solve = cloud is cloud_list[-1] # The last cloud computes the fields
                    spacech_ele.recompute_spchg_efield(MP_e, t_curr=beamtim.tt_curr, flag_solve=flag_solve, flag_reset=flag_reset)

                    # Copy rho to cloud
                    cloud.rho = spacech_ele.rho - sum([cl.rho for cl in cloud_list[:i_cloud]])


            for cloud in cloud_list:
                ## savings
                cloud.impact_man = cloud.pyeclsaver.witness(cloud.MP_e, beamtim, spacech_ele, cloud.impact_man, cloud.dynamics,
                                                            cloud.gas_ion_flag, cloud.resgasion, cloud.t_ion, t_sc_ON,
                                                            cloud.photoem_flag, cloud.phemiss, flag_presence_sec_beams,
                                                            sec_beams_list, cloud_list, rho_cloud = cloud.rho)

                ## every bunch passage
                if beamtim.flag_new_bunch_pass:

                    ## clean
                    cloud.MP_e.clean_small_MPs()

                    ## regeneration
                    cloud.MP_e.check_for_regeneration()

                    ## soft regeneration
                    cloud.MP_e.check_for_soft_regeneration()

            if beamtim.flag_new_bunch_pass:
                print '**** Done pass_numb = %d/%d\n'%(beamtim.pass_numb,beamtim.N_pass_tot)

            ## every bunch passage
            if t_end_sim is not None:
                if beamtim.tt_curr>    t_end_sim:
                    print 'Reached user defined t_end_sim --> Ending simulation'
                    break

    def load_state(self, filename_simulation_state, force_disable_save_simulation_state=True, filen_main_outp='Pyecltest_restarted') #, reset_pyeclsaver = True):

        print 'Realoading state from file: %s...'% filename_simulation_state

        with open(filename_simulation_state, 'rb') as fid:
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

            filen_outp_ext = filen_main_outp.split('Pyecltest')[-1]
            filen_outp_root = cloud.pyeclsaver.filen_main_outp.split('.mat')[0]
            cloud.pyeclsaver.filen_main_outp = filen_outp_root + filen_outp_ext + '.mat'

        print 'Restoring PyPIC LU object...'
        self.spacech_ele.PyPICobj.build_sparse_solver()
        print 'Done reload.'

#-Begin-preamble-------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 7.5.0
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

from __future__ import print_function
import scipy.io as sio
import numpy as np
import os
import subprocess
import hist_for as histf
import time
from scipy.constants import e as qe
import myfilemanager as mfm
try:
    # cPickle is faster in python2
    import cPickle as pickle
except ImportError:
    # No cPickle in python3
    import pickle


class pyecloud_saver:

    def __init__(self, logfile_path):
        print('Starting pyecloud_saver init.')
        self.logfile_path = logfile_path
        timestr = time.strftime("%d %b %Y %H:%M:%S", time.localtime())

        # These git commands return the hash and the branch of the specified git directory.
        path_to_git = os.path.dirname(os.path.abspath(__file__)) +'/.git'
        cmd_hash = 'git --git-dir %s rev-parse HEAD' % path_to_git
        cmd_branch = 'git --git-dir %s rev-parse --abbrev-ref HEAD' % path_to_git

        try:
            git_hash = 'git hash: %s' % (subprocess.check_output(cmd_hash.split()).split()[0])
        except Exception as e:
            git_hash = 'Retrieving git hash failed'
            print(e)
        print(git_hash)

        try:
            git_branch = 'git branch: %s' % (subprocess.check_output(cmd_branch.split()).split()[0])
        except Exception as e:
            git_branch = 'Retrieving git branch failed'
            print(e)
        print(git_branch)

        with open(self.logfile_path,'w') as flog:
            flog.write('PyECLOUD Version 7.5.0\n')
            flog.write('%s\n' % git_hash)
            flog.write('%s\n' % git_branch)
            flog.write('Simulation started on %s\n' % timestr)

        self.extract_sey = True


    def start_observing(self, Dt_ref, MP_e, beamtim, impact_man,
                 r_center, Dt_En_hist, logfile_path, progress_path, flag_detailed_MP_info =0,
                 cos_angle_width=0.05, flag_cos_angle_hist=True,
                 flag_movie=0, flag_sc_movie=0, save_mp_state_time_file=-1,
                 flag_presence_sec_beams=False, sec_beams_list=[], dec_fac_secbeam_prof=1,
                 el_density_probes = [],
                 save_simulation_state_time_file = -1,
                 x_min_hist_det=None, x_max_hist_det=None, y_min_hist_det=None, y_max_hist_det=None, Dx_hist_det=None,
                 filen_main_outp = 'Pyecltest', dec_fact_out = 1, stopfile = 'stop',
                 flag_multiple_clouds = False, cloud_name = None, flag_last_cloud = True,
                 checkpoint_DT=None):
        print('Start pyecloud_saver observation')

        self.filen_main_outp = filen_main_outp

        if '/' in self.filen_main_outp:
            self.folder_outp = '/'.join(self.filen_main_outp.split('/')[:-1])
        else:
            self.folder_outp='./'

        self.flag_detailed_MP_info = flag_detailed_MP_info

        # cloud info
        self.flag_multiple_clouds = flag_multiple_clouds
        self.cloud_name = cloud_name
        self.flag_last_cloud = flag_last_cloud

        # Init MP state saving
        self._MP_state_init(save_mp_state_time_file)

        # Init simulation state saving
        self._sim_state_init(save_simulation_state_time_file)

        # Init checkpoint saving
        self._checkpoint_init(checkpoint_DT)

        # Init charge distribution video saving
        self._rho_video_init(flag_movie)

        # Init electric field video saving
        self._sc_video_init(flag_sc_movie)

        # Init step by step data saving
        self._stepbystep_data_init(Dt_ref, dec_fact_out, el_density_probes, r_center,
                                    initial_size_t_vect=1000)

        # Init pass by pass data saving
        self._pass_by_pass_data_init(impact_man,
                    x_min_hist_det, x_max_hist_det, y_min_hist_det, y_max_hist_det,  Dx_hist_det)

        # Init energy and cos angle histogram saving
        self._energy_and_cos_angle_hist_init(Dt_En_hist, flag_cos_angle_hist, cos_angle_width)

        #Space charge electrostatic energy
        self.t_sc_video=[]
        self.U_sc_eV=[]

        # Prepare space for density histogram calculation
        self.nel_hist_line = np.zeros(impact_man.Nxg_hist,float)

        #logfile and progress file
        self.logfile_path=logfile_path
        self.progress_path=progress_path
        self.stopfile = stopfile

        # Store secondary beam profiles
        if flag_presence_sec_beams:
            self.t_sec_beams = beamtim.t[::dec_fac_secbeam_prof]
            self.sec_beam_profiles = []
            for sec_beam in sec_beams_list:
                self.sec_beam_profiles.append(sec_beam.lam_t_array[::dec_fac_secbeam_prof])
            self.sec_beam_profiles=np.array(self.sec_beam_profiles)
        else:
            self.t_sec_beams = -1
            self.sec_beam_profiles = -1


        # extract SEY curves
        if self.extract_sey:
            n_rep = 10000
            self.sey_test_E_impact_eV = np.array(list(np.arange(0, 499., 1.)) + list(np.arange(500., 2000, 5)))
            self.sey_test_cos_theta = np.linspace(0, 1., 10)
            self.sey_test_del_true_mat, self.sey_test_del_elast_mat = impact_man.extract_sey_curves(n_rep, self.sey_test_E_impact_eV, self.sey_test_cos_theta)
        else:
            self.sey_test_E_impact_eV = 0.
            self.sey_test_cos_theta = 0.
            self.sey_test_del_true_mat, self.sey_test_del_elast_mat = 0.,0.


        if checkpoint_DT is not None:
            # blaaaaa
            pass

        # Log
        print('Done init pyecloud_saver.')
        flog=open(self.logfile_path,'a')
        timestr = time.strftime("%d %b %Y %H:%M:%S", time.localtime())
        flog.write('Initialization finished on %s\n'%timestr)
        flog.close()

    def witness(self, MP_e, beamtim, spacech_ele, impact_man,
                dynamics,gas_ion_flag,resgasion,t_ion,
                t_sc_ON, photoem_flag, phemiss,flag_presence_sec_beams,sec_beams_list,
                cloud_list, rho_cloud = None):

        ####################################################
        # Quantites saved at custom times provided by user #
        ####################################################

        # Check for MP save state
        self._MP_state_save(MP_e, beamtim)

        # Check for simulation save state
        self._sim_state_save(beamtim, spacech_ele, t_sc_ON, flag_presence_sec_beams,
                              sec_beams_list, self.flag_multiple_clouds, cloud_list)

        # Check for checkpoint save state
        self._checkpoint_save(beamtim, spacech_ele, t_sc_ON, flag_presence_sec_beams,
                            sec_beams_list, self.flag_multiple_clouds, cloud_list)

        # Check for save video charge density
        self._rho_video_save(spacech_ele, beamtim, rho_cloud)

        # Check for save video electric field
        self._sc_video_save(spacech_ele, beamtim)

        # Check for energy and cos angle hist update
        self._energy_and_cos_angle_hist_save(beamtim, impact_man)

        #Space charge electrostatic energy
        if spacech_ele.last_recomputation_check:
            self.t_sc_video.append(beamtim.tt_curr)
            self.U_sc_eV.append(spacech_ele.U_sc_eV_stp)

        #########################
        # Step-by-step saveings #
        #########################
        self._stepbystep_data_save(impact_man, MP_e, beamtim)


        ##########################################################
        # Quantites saved at each bunch passage and dump to file #
        ##########################################################

        if beamtim.flag_new_bunch_pass:

            self._pass_by_pass_data_save(MP_e, impact_man, beamtim)

            # I want all elements that go to the output file to be members of this object
            self.xg_hist = impact_man.xg_hist
            self.En_g_hist = impact_man.En_g_hist
            self.b_spac = beamtim.b_spac
            self.area = impact_man.chamb.area

            sio.savemat(self.filen_main_outp, self.build_outp_dict(), oned_as='row')

        if beamtim.flag_new_bunch_pass:
            self._logfile_progressfile_stofile(beamtim, MP_e)


        return impact_man

    def _pass_by_pass_data_init(self, impact_man,
                    x_min_hist_det, x_max_hist_det, y_min_hist_det, y_max_hist_det,  Dx_hist_det):

        #pass by pass data
        self.t_hist=[]
        self.nel_impact_hist_tot = []
        self.nel_impact_hist_scrub = []
        self.energ_eV_impact_hist = []
        self.nel_hist = []
        self.N_mp_impact_pass = []
        self.N_mp_corrected_pass = []
        self.N_mp_pass = []
        self.N_mp_ref_pass = []

        if impact_man.flag_seg:
                self.nel_hist_impact_seg = []
                self.energ_eV_impact_seg = []
        else:
                self.nel_hist_impact_seg = -1
                self.energ_eV_impact_seg = -1

        # detailed hist
        self.flag_hist_det = False
        self.xg_hist_det = -1
        self.nel_hist_det = -1
        if x_min_hist_det is not None:
            if x_max_hist_det is None or y_min_hist_det is None or  y_max_hist_det is None or  Dx_hist_det is None:
                raise ValueError('x_min_hist_det is not None but one among x_max_hist_det, y_min_hist_det, y_max_hist_det, Dx_hist_det is None!')

            self.flag_hist_det = True
            self.x_min_hist_det = x_min_hist_det
            self.x_max_hist_det = x_max_hist_det
            self.y_min_hist_det = y_min_hist_det
            self.y_max_hist_det = y_max_hist_det
            self.Dx_hist_det = Dx_hist_det

            self.xg_hist_det=np.arange(x_min_hist_det-Dx_hist_det,x_max_hist_det+Dx_hist_det+0.1*Dx_hist_det,Dx_hist_det,float)
            self.Nxg_hist_det=len(self.xg_hist_det);
            self.bias_x_hist_det=min(self.xg_hist_det);

            self.nel_hist_det_line=np.zeros(self.Nxg_hist_det,float)
            self.nel_hist_det=[]

    def _pass_by_pass_data_save(self, MP_e, impact_man, beamtim):
        #update histograms
        self.nel_hist_line=0.0*self.nel_hist_line
        if MP_e.N_mp>0:
            histf.compute_hist(MP_e.x_mp[0:MP_e.N_mp],MP_e.nel_mp[0:MP_e.N_mp],impact_man.bias_x_hist,impact_man.Dx_hist,self.nel_hist_line)


            # detailed histogram
            if self.flag_hist_det:
                #print 'here 1'
                mask_det_hist_x = np.logical_and(MP_e.x_mp[0:MP_e.N_mp]>self.x_min_hist_det, MP_e.x_mp[0:MP_e.N_mp]<self.x_max_hist_det)
                mask_det_hist_y = np.logical_and(MP_e.y_mp[0:MP_e.N_mp]>self.y_min_hist_det, MP_e.y_mp[0:MP_e.N_mp]<self.y_max_hist_det)
                mask_det_hist = np.logical_and(mask_det_hist_x, mask_det_hist_y)

                self.nel_hist_det_line=0.0*self.nel_hist_det_line

                if np.sum(mask_det_hist)>0:
                    #print 'here 2'
                    histf.compute_hist(MP_e.x_mp[0:MP_e.N_mp][mask_det_hist],MP_e.nel_mp[0:MP_e.N_mp][mask_det_hist],
                                       self.bias_x_hist_det, self.Dx_hist_det, self.nel_hist_det_line)


        self.nel_hist.append(self.nel_hist_line.copy())
        self.t_hist.append(beamtim.tt_curr)
        self.nel_impact_hist_tot.append(impact_man.nel_impact_hist_tot.copy())
        impact_man.reset_impact_hist_tot()
        self.nel_impact_hist_scrub.append(impact_man.nel_impact_hist_scrub.copy())
        impact_man.reset_impact_hist_scrub()
        self.energ_eV_impact_hist.append(impact_man.energ_eV_impact_hist.copy())
        impact_man.reset_energ_eV_impact_hist()

        self.N_mp_impact_pass.append(impact_man.chamb.N_mp_impact)
        self.N_mp_corrected_pass.append(impact_man.chamb.N_mp_corrected)
        self.N_mp_pass.append(MP_e.N_mp)
        self.N_mp_ref_pass.append(MP_e.nel_mp_ref)

        if impact_man.flag_seg:
            self.nel_hist_impact_seg.append(impact_man.nel_hist_impact_seg.copy())
            impact_man.reset_hist_impact_seg()

        if impact_man.flag_seg:
            self.energ_eV_impact_seg.append(impact_man.energ_eV_impact_seg.copy())
            impact_man.reset_energ_impact_seg()

        if self.flag_hist_det:
            self.nel_hist_det.append(self.nel_hist_det_line.copy())

    def build_outp_dict(self):
        saved_dict = {
                    't_hist':               self.t_hist,
                    'nel_hist':             self.nel_hist,
                    'xg_hist':              self.xg_hist,
                    'nel_impact_hist_tot':  self.nel_impact_hist_tot,
                    'nel_impact_hist_scrub':self.nel_impact_hist_scrub,
                    'energ_eV_impact_hist': self.energ_eV_impact_hist,
                    'En_g_hist':            self.En_g_hist,
                    'En_hist':              self.En_hist,
                    't_En_hist':            self.t_En_hist,
                    'b_spac':               self.b_spac,
                    't_sc_video':           np.array(self.t_sc_video),
                    'U_sc_eV':              np.array(self.U_sc_eV),
                    'N_mp_impact_pass':     self.N_mp_impact_pass,
                    'N_mp_corrected_pass':  self.N_mp_corrected_pass,
                    'N_mp_pass':            self.N_mp_pass,
                    'N_mp_ref_pass':        self.N_mp_ref_pass,
                    'nel_hist_impact_seg':  self.nel_hist_impact_seg,
                    'energ_eV_impact_seg':  self.energ_eV_impact_seg,
                    't_sec_beams':          self.t_sec_beams,
                    'sec_beam_profiles':    self.sec_beam_profiles,
                    'x_el_dens_probes':     self.x_el_dens_probes,
                    'y_el_dens_probes':     self.y_el_dens_probes,
                    'r_el_dens_probes':     self.r_el_dens_probes,
                    'nel_hist_det':         self.nel_hist_det,
                    'xg_hist_det':          self.xg_hist_det,
                    'dec_fact_out':         self.dec_fact_out,
                    'sey_test_E_impact_eV': self.sey_test_E_impact_eV,
                    'sey_test_cos_theta':   self.sey_test_cos_theta,
                    'sey_test_del_true_mat':self.sey_test_del_true_mat,
                    'sey_test_del_elast_mat':self.sey_test_del_elast_mat,
                    'chamber_area':         self.area,
                    'cos_angle_hist':       self.cos_angle_hist,
                    'xg_hist_cos_angle':    self.xg_hist_cos_angle
                }

        saved_dict.update(self._stepbystep_get_dict())

        for kk in saved_dict.keys():
            saved_dict[kk] = np.array(saved_dict[kk])

        return saved_dict

    def _checkpoint_init(self, checkpoint_DT):
        # Simulation state saver init
        if checkpoint_DT is None:
            self.flag_save_checkpoint=False
        elif type(checkpoint_DT) is int and checkpoint_DT == -1:
            self.flag_save_checkpoint=False
        else:
            self.flag_save_checkpoint=True
            self.checkpoint_DT = checkpoint_DT
            self.t_last_checkp = 0
            self.i_checkp = 0

    def _checkpoint_save(self, beamtim, spacech_ele, t_sc_ON, flag_presence_sec_beams,
                    sec_beams_list, flag_multiple_clouds, cloud_list):
        # First check if it is time to save a checkpoint
        if (self.flag_save_checkpoint):
            if (beamtim.tt_curr - self.t_last_checkp >= self.checkpoint_DT):
                
                outpath = self.folder_outp + 'simulation_checkpoint_%d.pkl'%(self.i_checkp)
            
                self._sim_state_single_save(beamtim, spacech_ele, t_sc_ON, flag_presence_sec_beams,
                            sec_beams_list, flag_multiple_clouds, cloud_list, outpath)
                print('Save simulation checkpoint in: ' + outpath)

                if self.i_checkp>0:
                    prevpath = outpath = self.folder_outp + 'simulation_checkpoint_%d.pkl'%(self.i_checkp-1)
                    os.remove(prevpath)
                    print('Removed simulation checkpoint in: ' + prevpath)                        
                        
                self.i_checkp += 1
                self.t_last_checkp = beamtim.tt_curr
                

    def load_from_output(self, fname, last_t=None):

        #restore the Pyecltest.mat up to last t
        ob = mfm.myloadmat_to_obj(self.folder_outp + '/' + fname)
        dict_history = mfm.obj_to_dict(ob)

        idx_t = (np.abs(dict_history['t'] - last_t)).argmin()  # index closest to last_t

        idx_t_hist = (np.abs(dict_history['t_hist'] - last_t)).argmin()

        # Delete everything in Pyecltest.mat recorded after the last checkpoint
        saved_every_timestep_list = ['En_emit_eV_time',
                                     'En_imp_eV_time',
                                     'En_kin_eV_time',
                                     'Nel_emit_time',
                                     'Nel_imp_time',
                                     'Nel_timep',
                                     'cen_density',
                                     'lam_t_array',
                                     'N_mp_time',
                                     't']

        saved_every_passage_list = ['En_hist',
                                    't_En_hist',
                                    'N_mp_corrected_pass',
                                    'N_mp_impact_pass',
                                    'N_mp_pass',
                                    'N_mp_ref_pass',
                                    'energ_eV_impact_hist',
                                    'energ_eV_impact_seg',
                                    'nel_hist',
                                    'nel_hist_det',
                                    'nel_hist_impact_seg',
                                    'nel_impact_hist_scrub',
                                    'nel_impact_hist_tot',
                                    'cos_angle_hist',
                                    't_hist']

        not_time_dependent_list = ['xg_hist',
                                   'xg_hist_det',
                                   'En_g_hist',
                                   'b_spac',
                                   't_sc_video',
                                   't_sec_beams',
                                   'sec_beam_profiles',
                                   'x_el_dens_probes',
                                   'y_el_dens_probes',
                                   'r_el_dens_probes',
                                   'xg_hist_cos_angle',
                                   'el_dens_at_probes',
                                   'dec_fact_out',
                                   'chamber_area',
                                   'sey_test_del_true_mat',
                                   'sey_test_E_impact_eV',
                                   'sey_test_del_elast_mat',
                                   'sey_test_cos_theta',
                                   'U_sc_eV']

        should_be_list_list = ['t_sc_video',
                               'U_sc_eV',
                               'x_el_dens_probes',
                               'y_el_dens_probes',
                               'r_el_dens_probes',
                                ]
        dict_restored = {}
        for var in saved_every_timestep_list:
            if var in dict_history.keys():
                if dict_history[var].shape == np.array(0).shape:
                    dict_restored[var] = dict_history[var]
                else:
                    dict_restored[var] = dict_history[var][: idx_t]

        for var in saved_every_passage_list:
            if var in dict_history.keys():
                if dict_history[var].shape == np.array(0).shape:
                    dict_restored[var] = dict_history[var].tolist()
                else:
                    dict_restored[var] = dict_history[var][: idx_t_hist].tolist()

        for var in not_time_dependent_list:
            if var in dict_history.keys():
                if var in should_be_list_list:
                    dict_restored[var] = dict_history[var].tolist()
                else:
                    dict_restored[var] = dict_history[var]

        # Restore this pyecloud_saver object with values from dict_restored
        for var in dict_restored.keys():
            setattr(self, var, dict_restored[var])


    def _stepbystep_check_for_data_resize(self):
        if self.i_last_save==(len(self.t_dec)-1):
            print('Saver: resizing from %d to %d...'%(len(self.t_dec), 2*len(self.t_dec)))
            list_members = [
                't_dec',
                'lam_t_array_dec',
                'Nel_time',
                'Nel_imp_time',
                'Nel_emit_time',
                'En_imp_eV_time',
                'En_emit_eV_time',
                'En_kin_eV_time',
                'cen_density'
                ]
            if self.flag_detailed_MP_info==1:
                list_members.append('N_mp_time')

            for mm in list_members:
                vv = getattr(self, mm)
                setattr(self, mm, np.concatenate((vv, 0*vv)))

            if self.flag_el_dens_probes:
                self.el_dens_at_probes = np.concatenate(
                    (self.el_dens_at_probes, 0*self.el_dens_at_probes), axis=1)
            print('Done resizing')

    def _stepbystep_data_init(self, Dt_ref, dec_fact_out, el_density_probes, r_center, initial_size_t_vect):

        #step by step data

        # introduce decimation
        self.Dt_ref = Dt_ref
        self.dec_fact_out = dec_fact_out
        self.Dt_save = (dec_fact_out-0.0001)*Dt_ref
        self.i_last_save = -1
        self.t_last_save = -1.

        self.r_center=r_center

        self.Nel_impact_last_step_group = 0
        self.Nel_emit_last_step_group = 0
        self.En_imp_last_step_group_eV = 0
        self.En_emit_last_step_group_eV = 0

        #step by step data
        self.t_dec = np.zeros(initial_size_t_vect, dtype=float)
        self.lam_t_array_dec = 0*self.t_dec
        self.Nel_time=0.*self.t_dec;
        self.Nel_imp_time=0.*self.t_dec;
        self.Nel_emit_time=0.*self.t_dec;
        self.En_imp_eV_time=0.*self.t_dec;
        self.En_emit_eV_time=0.*self.t_dec;
        self.En_kin_eV_time=0.*self.t_dec;
        self.cen_density=0.*self.t_dec;


        if self.flag_detailed_MP_info==1:
            self.N_mp_time=0.*self.t_dec
        else:
            self.N_mp_time=-1

        # initialize electron density probes
        self.flag_el_dens_probes = False
        self.x_el_dens_probes = -1
        self.y_el_dens_probes = -1
        self.r_el_dens_probes = -1
        self.el_dens_at_probes = -1
        if len(el_density_probes)>0:
            self.flag_el_dens_probes = True
            self.N_el_dens_probes = len(el_density_probes)
            self.el_dens_at_probes = np.zeros((self.N_el_dens_probes, len(self.t_dec))) # to be changed
            self.x_el_dens_probes = []
            self.y_el_dens_probes = []
            self.r_el_dens_probes = []
            for ii in xrange(self.N_el_dens_probes):
                self.x_el_dens_probes.append(el_density_probes[ii]['x'])
                self.y_el_dens_probes.append(el_density_probes[ii]['y'])
                self.r_el_dens_probes.append(el_density_probes[ii]['r_obs'])

            self.x_el_dens_probes = np.array(self.x_el_dens_probes)
            self.y_el_dens_probes = np.array(self.y_el_dens_probes)
            self.r_el_dens_probes = np.array(self.r_el_dens_probes)



    def _stepbystep_data_save(self, impact_man, MP_e, beamtim):
        #save step by step data
        # Vars to be accumulated
        self.Nel_impact_last_step_group += impact_man.Nel_impact_last_step
        self.Nel_emit_last_step_group += impact_man.Nel_emit_last_step
        self.En_imp_last_step_group_eV += impact_man.En_imp_last_step_eV
        self.En_emit_last_step_group_eV += impact_man.En_emit_last_step_eV


        #if np.mod(beamtim.ii_curr, self.dec_fact_out)==0:
        if beamtim.tt_curr-self.t_last_save >= self.Dt_save:


            self._stepbystep_check_for_data_resize()

            self.i_last_save+=1
            self.t_last_save = beamtim.tt_curr

            self.t_dec[self.i_last_save] = beamtim.tt_curr
            self.lam_t_array_dec[self.i_last_save] = beamtim.lam_t_curr

            self.Nel_imp_time[self.i_last_save] = self.Nel_impact_last_step_group
            self.Nel_emit_time[self.i_last_save] = self.Nel_emit_last_step_group
            self.En_imp_eV_time[self.i_last_save] = self.En_imp_last_step_group_eV
            self.En_emit_eV_time[self.i_last_save] = self.En_emit_last_step_group_eV

            self.Nel_impact_last_step_group = 0
            self.Nel_emit_last_step_group = 0
            self.En_imp_last_step_group_eV = 0
            self.En_emit_last_step_group_eV = 0


            self.Nel_time[self.i_last_save]=np.sum(MP_e.nel_mp[0:MP_e.N_mp]);
            self.En_kin_eV_time[self.i_last_save]=np.sum(0.5*MP_e.mass/qe*MP_e.nel_mp[0:MP_e.N_mp]*(MP_e.vx_mp[0:MP_e.N_mp]*MP_e.vx_mp[0:MP_e.N_mp]+MP_e.vy_mp[0:MP_e.N_mp]*MP_e.vy_mp[0:MP_e.N_mp]+MP_e.vz_mp[0:MP_e.N_mp]*MP_e.vz_mp[0:MP_e.N_mp]));

            flag_center=((MP_e.x_mp**2 + MP_e.y_mp**2)<self.r_center**2);
            flag_center[MP_e.N_mp:]=False
            self.cen_density[self.i_last_save]=np.sum(MP_e.nel_mp[flag_center])/(np.pi*self.r_center*self.r_center)

            if self.flag_el_dens_probes:
                for ii in xrange(self.N_el_dens_probes):
                    flag_center=((MP_e.x_mp-self.x_el_dens_probes[ii])**2 + (MP_e.y_mp-self.y_el_dens_probes[ii])**2)<self.r_el_dens_probes[ii]**2;
                    flag_center[MP_e.N_mp:]=False
                    self.el_dens_at_probes[ii, self.i_last_save]=np.sum(MP_e.nel_mp[flag_center])/(np.pi*self.r_el_dens_probes[ii]**2)


            if self.flag_detailed_MP_info==1:
                self.N_mp_time[self.i_last_save]=MP_e.N_mp

    def _stepbystep_get_dict(self):

        dict_sbs_data = {
            't':self.t_dec[:self.i_last_save+1],
            'lam_t_array':self.lam_t_array_dec[:self.i_last_save+1],
            'Nel_timep':self.Nel_time[:self.i_last_save+1],
            'Nel_imp_time':self.Nel_imp_time[:self.i_last_save+1],
            'Nel_emit_time':self.Nel_emit_time[:self.i_last_save+1],
            'En_imp_eV_time':self.En_imp_eV_time[:self.i_last_save+1],
            'En_emit_eV_time':self.En_emit_eV_time[:self.i_last_save+1],
            'En_kin_eV_time':self.En_kin_eV_time[:self.i_last_save+1],
            'cen_density':self.cen_density[:self.i_last_save+1]
        }

        if self.flag_detailed_MP_info==1:
            dict_sbs_data['N_mp_time'] = self.N_mp_time[:self.i_last_save+1]

        if self.flag_el_dens_probes:
            dict_sbs_data['el_dens_at_probes'] = self.el_dens_at_probes[:, :self.i_last_save]

        return dict_sbs_data


    def _MP_state_init(self, save_mp_state_time_file):
        # MP state saver init
        try:
            save_mp_state_time_file[0] #check if iterable
            self.flag_save_MP_state=True
            if type(save_mp_state_time_file) is str:
                dict_save_mp_state_time=sio.loadmat(save_mp_state_time_file)
                self.t_obs=np.squeeze(dict_save_mp_state_time['t_obs'].real)
            else:
                self.t_obs = np.array(save_mp_state_time_file)

            self.N_obs=len(self.t_obs)
            self.i_obs=0
        except TypeError:
            self.flag_save_MP_state=False

    def _MP_state_save(self, MP_e, beamtim):
        #MP state save
        if self.flag_save_MP_state:
            if  (MP_e.N_mp>0) and (self.i_obs<self.N_obs):
                if (beamtim.tt_curr>=self.t_obs[self.i_obs]):
                    if self.flag_multiple_clouds:
                        filename_MP_state='MP_state_%s_%d'%(self.cloud_name,self.i_obs)
                    else:
                        filename_MP_state='MP_state_%d'%(self.i_obs)
                    path_MP_state = self.folder_outp+'/'+filename_MP_state
                    sio.savemat(path_MP_state,{'tt':beamtim.tt_curr,'N_mp':MP_e.N_mp, 'x_mp':MP_e.x_mp[0:MP_e.N_mp], 'y_mp':MP_e.y_mp[0:MP_e.N_mp], 'z_mp':MP_e.z_mp[0:MP_e.N_mp],\
                                                    'vx_mp':MP_e.vx_mp[0:MP_e.N_mp], 'vy_mp':MP_e.vy_mp[0:MP_e.N_mp], 'vz_mp':MP_e.vz_mp[0:MP_e.N_mp], 'nel_mp':MP_e.nel_mp[0:MP_e.N_mp]},oned_as='row')

                    print('Save MP state in: ' + path_MP_state)
                    self.i_obs=self.i_obs+1


    def _sim_state_init(self, save_simulation_state_time_file):
        # Simulation state saver init
        if save_simulation_state_time_file is None:
            self.flag_save_simulation_state=False
        elif type(save_simulation_state_time_file) is int and save_simulation_state_time_file == -1:
            self.flag_save_simulation_state=False
        else:
            self.flag_save_simulation_state=True
            if type(save_simulation_state_time_file) is str:
                dict_save_simulation_state_time=sio.loadmat(save_simulation_state_time_file)
                self.t_obs_sim=np.squeeze(dict_save_simulation_state_time['t_obs'].real)
            else:
                self.t_obs_sim = np.array(save_simulation_state_time_file)

            self.N_obs_sim=len(self.t_obs_sim)
            self.i_obs_sim=0
            
    def _sim_state_single_save(self, beamtim, spacech_ele, t_sc_ON, flag_presence_sec_beams,
                    sec_beams_list, flag_multiple_clouds, cloud_list, outfile):
        
        if self.flag_last_cloud:
            
            temp_luobj = spacech_ele.PyPICobj.luobj
            spacech_ele.luobj=None
            spacech_ele.PyPICobj.luobj=None

            #~ dynamics.get_B=None

            # remove savers
            temp_saver_list = []
            for cloud in cloud_list:
                temp_saver_list.append(cloud.pyeclsaver)
                cloud.pyeclsaver = 'removed'

            dict_state = {
            'beamtim':beamtim,
            'spacech_ele':spacech_ele,
            't_sc_ON':t_sc_ON,
            'flag_presence_sec_beams':flag_presence_sec_beams,
            'sec_beams_list':sec_beams_list,
            'flag_multiple_clouds':self.flag_multiple_clouds,
            'cloud_list':cloud_list}

            with open(outfile, 'wb') as fid:
                # use best protocol available
                pickle.dump(dict_state, fid, protocol=-1)

            # put back savers
            for cloud, saver in zip(cloud_list, temp_saver_list):
                cloud.pyeclsaver = saver

            spacech_ele.PyPICobj.luobj = temp_luobj

            print('Save simulation state in: ' + outfile)

    def _sim_state_save(self, beamtim, spacech_ele, t_sc_ON, flag_presence_sec_beams,
                    sec_beams_list, flag_multiple_clouds, cloud_list):
        #Simulation state save
        if self.flag_save_simulation_state:
            if self.i_obs_sim<self.N_obs_sim:
                if (beamtim.tt_curr>=self.t_obs_sim[self.i_obs_sim]):
                    filename_simulation_state='simulation_state_%d.pkl'%(self.i_obs_sim)
                    outpath = self.folder_outp+'/'+filename_simulation_state
                    
                    self._sim_state_single_save(beamtim, spacech_ele, t_sc_ON, flag_presence_sec_beams,
                        sec_beams_list, flag_multiple_clouds, cloud_list, outpath)
                    
                    self.i_obs_sim=self.i_obs_sim+1

    def _rho_video_init(self, flag_movie):
        #rho video
        self.flag_video=(flag_movie==1)
        self.rho_video=None
        self.t_video=None
        #rho video cloud
        self.rho_video_cloud=None
        self.t_video_cloud=None

    def _rho_video_save(self, spacech_ele, beamtim, rho_cloud):
        #save rho video
        if self.flag_video and self.flag_last_cloud:
            if not os.path.exists(self.folder_outp+'/rho_video'):
                os.makedirs(self.folder_outp+'/rho_video')
            if self.rho_video is None:
                self.rho_video=[]
                self.t_video=[]
            if spacech_ele.last_recomputation_check:
                self.rho_video.append(spacech_ele.rho)
                self.t_video.append(beamtim.tt_curr)
            if beamtim.flag_new_bunch_pass:
                self.rho_video=np.array(self.rho_video)
                self.t_video=np.array(self.t_video)
                filename_rho=self.folder_outp+'/rho_video/rho_pass%d.mat'%(beamtim.pass_numb-1)
                print('Saving %s'%filename_rho)
                sio.savemat(filename_rho,{'xg_sc':spacech_ele.xg,'yg_sc':spacech_ele.yg,'t_video':self.t_video,'rho_video':self.rho_video},oned_as='row')
                print('Done')
                self.rho_video=[]
                self.t_video=[]

        # save rho video for cloud
        if self.flag_video and self.flag_multiple_clouds:
            if not os.path.exists(self.folder_outp+'/rho_video_%s'%(self.cloud_name)):
                os.makedirs(self.folder_outp+'rho_video_%s'%(self.cloud_name))
            if self.rho_video_cloud is None:
                    self.rho_video_cloud = []
                    self.t_video_cloud = []
            if spacech_ele.last_recomputation_check:
                if rho_cloud is None:
                    print('Warning! No rho provided for saving.')
                else:
                    self.rho_video_cloud.append(rho_cloud)
                self.t_video_cloud.append(beamtim.tt_curr)
            if beamtim.flag_new_bunch_pass:
                self.rho_video_cloud=np.array(self.rho_video_cloud)
                self.t_video_cloud=np.array(self.t_video_cloud)
                filename_rho=self.folder_outp+'rho_video_%s/rho_pass%d.mat'%(self.cloud_name, beamtim.pass_numb-1)
                print('Saving %s'%filename_rho)
                sio.savemat(filename_rho,{'xg_sc':spacech_ele.xg,'yg_sc':spacech_ele.yg,'t_video':self.t_video_cloud,'rho_video':self.rho_video_cloud},oned_as='row')
                print('Done')
                self.rho_video_cloud=[]
                self.t_video_cloud=[]

    def _sc_video_init(self, flag_sc_movie):
        #efield video
        self.flag_sc_video=(flag_sc_movie==1)
        self.efx_video=None
        self.efy_video=None
        self.t_efield_video=None

    def _sc_video_save(self, spacech_ele, beamtim):
        #save efield video
        if self.flag_sc_video:
            if not os.path.exists(self.folder_outp+'/efield_video'):
                os.makedirs(self.folder_outp+'/efield_video')
            if self.efx_video is None:
                self.efx_video=[]
                self.efy_video=[]
                self.t_efield_video=[]
            if spacech_ele.last_recomputation_check:
                self.efx_video.append(spacech_ele.efx)
                self.efy_video.append(spacech_ele.efy)
                self.t_efield_video.append(beamtim.tt_curr)
            if beamtim.flag_new_bunch_pass:
                self.efx_video=np.array(self.efx_video)
                self.efy_video=np.array(self.efy_video)
                self.t_efield_video=np.array(self.t_efield_video)
                filename_efield=self.folder_outp+'/efield_video/efield_pass%d.mat'%(beamtim.pass_numb-1)
                print('Saving %s'%filename_efield)
                sio.savemat(filename_efield,{'xg_sc':spacech_ele.xg,'yg_sc':spacech_ele.yg,'t_efield_video':self.t_efield_video,
                                          'efx_video':self.efx_video, 'efy_video':self.efy_video},oned_as='row')
                print('Done')
                self.efx_video=[]
                self.efy_video=[]
                self.t_efield_video=[]

    def _logfile_progressfile_stofile(self, beamtim, MP_e):
            # logfile and progressfile
            timestr = time.strftime("%d %b %Y %H:%M:%S", time.localtime())

            string_tolog= timestr+(' pass. %d/%d, cloud=%s: Nel_tot=%e N_mp=%d\n'%(beamtim.pass_numb,beamtim.N_pass_tot,self.cloud_name,np.sum(MP_e.nel_mp[0:MP_e.N_mp]),MP_e.N_mp))

            try:
                with open(self.logfile_path,'a') as flog:
                    flog.write(string_tolog)
            except IOError as err:
                print('Got: ',err)
                print('while trying to write the following line on logfile:')
                print(string_tolog)

            try:
                with open(self.progress_path,'w') as flog:
                    flog.write(('%f'%(float(beamtim.pass_numb)/float(beamtim.N_pass_tot))))
            except IOError as err:
                print('Got: ',err)
                print('while trying to write the following line on progress file:')
                print('%f'%(float(beamtim.pass_numb)/float(beamtim.N_pass_tot)))

            #stop simulation
            try:
                with open(self.stopfile, 'r') as fid:
                    fid.readlines()
                    raise ValueError('Stopped by user.')
            except IOError:
                pass

    def _energy_and_cos_angle_hist_init(self, Dt_En_hist, flag_cos_angle_hist,
                cos_angle_width):
        # Energy histogram init
        self.Dt_En_hist = Dt_En_hist
        self.t_last_En_hist = -1.
        self.En_hist=[]
        self.t_En_hist = []

        # Angle histogram
        self.flag_cos_angle_hist = flag_cos_angle_hist
        if flag_cos_angle_hist:
            N_angles = int(1./ cos_angle_width)+1
            self.cos_angle_hist = []
            self.xg_hist_cos_angle = np.linspace(0., 1., N_angles)
        else:
            self.cos_angle_hist = -1
            self.xg_hist_cos_angle = -1

    def _energy_and_cos_angle_hist_save(self, beamtim, impact_man):
        # Energy histogram saver
        # if (np.mod(beamtim.ii_curr,self.Nst_En_hist)==0):
        if beamtim.tt_curr>=self.t_last_En_hist+self.Dt_En_hist:
            self.En_hist.append(impact_man.En_hist_line.copy())
            self.t_En_hist.append(beamtim.tt_curr)
            impact_man.reset_En_hist_line()
            self.t_last_En_hist = beamtim.tt_curr

            if self.flag_cos_angle_hist:
                self.cos_angle_hist.append(impact_man.cos_angle_hist.copy())
                impact_man.reset_cos_angle_hist()

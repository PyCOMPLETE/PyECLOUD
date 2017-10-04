#----------------------------------------------------------------------
#----------------------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 6.4.1
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

from __future__ import print_function
import scipy.io as sio
import numpy as np
import os
import subprocess
import hist_for as histf
import time
try:
    # cPickle is faster in python2
    import cPickle as pickle
except ImportError:
    # No cPickle in python3
    import pickle


me=9.10938291e-31
qe=1.602176565e-19
qm=qe/me

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
            flog.write('PyECLOUD Version 6.4.1\n')
            flog.write('%s\n' % git_hash)
            flog.write('%s\n' % git_branch)
            flog.write('Simulation started on %s\n' % timestr)


    def start_observing(self, MP_e, beamtim, impact_man,
                 r_center, Dt_En_hist, logfile_path, progress_path, flag_detailed_MP_info=0,
                 cos_angle_width=0.05, flag_cos_angle_hist=True,
                 flag_movie=0, flag_sc_movie=0, save_mp_state_time_file=-1,
                 flag_presence_sec_beams=False, sec_beams_list=[], dec_fac_secbeam_prof=1,
                 el_density_probes = [],
                 save_simulation_state_time_file = -1,
                 x_min_hist_det=None, x_max_hist_det=None, y_min_hist_det=None, y_max_hist_det=None, Dx_hist_det=None,
                 filen_main_outp = 'Pyecltest', dec_fact_out = 1, stopfile = 'stop'):
        print('Start pyecloud_saver observation')

        self.filen_main_outp = filen_main_outp


        # introduce decimation
        self.dec_fact_out = dec_fact_out
        self.t_dec = beamtim.t[::dec_fact_out]
        self.lam_t_array_dec = beamtim.lam_t_array[::dec_fact_out]
        self.Nel_impact_last_step_group = 0
        self.Nel_emit_last_step_group = 0
        self.En_imp_last_step_group_eV = 0
        self.En_emit_last_step_group_eV = 0

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


        # Simulation state saver init
        if type(save_simulation_state_time_file) is int and save_simulation_state_time_file == -1:
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

        # Energy histogram init
        self.Nst_En_hist=int(round(Dt_En_hist/beamtim.Dt)) #number of steps per hist. line
        self.N_En_hist=int(float(beamtim.Nt)/float(self.Nst_En_hist))+2    #hist size in time dimension
        self.t_En_hist=np.zeros(self.N_En_hist,float)
        self.En_hist=np.zeros((self.N_En_hist,impact_man.Nbin_En_hist),float) #allocate histograms
        self.i_En_hist=0

        # Angle histogram
        self.flag_cos_angle_hist = flag_cos_angle_hist
        if flag_cos_angle_hist:
            N_angles = int(1./ cos_angle_width)+1
            self.cos_angle_hist = np.zeros((self.N_En_hist, N_angles),float)
            self.xg_hist_cos_angle = np.linspace(0., 1., N_angles)

        #Space charge electrostatic energy
        self.t_sc_video=[]
        self.U_sc_eV=[]

        #rho video
        self.flag_video=(flag_movie==1)
        self.rho_video=None
        self.t_video=None

        #efield video
        self.flag_sc_video=(flag_sc_movie==1)
        self.efx_video=None
        self.efy_video=None
        self.t_efield_video=None

        #step by step data
        self.r_center=r_center
        self.Nel_time=0.*self.t_dec;
        self.Nel_imp_time=0.*self.t_dec;
        self.Nel_emit_time=0.*self.t_dec;
        self.En_imp_eV_time=0.*self.t_dec;
        self.En_emit_eV_time=0.*self.t_dec;
        self.En_kin_eV_time=0.*self.t_dec;
        self.cen_density=0.*self.t_dec;

        self.flag_detailed_MP_info=flag_detailed_MP_info

        if flag_detailed_MP_info==1:
            self.N_mp_time=0.*self.t_dec
        else:
            self.N_mp_time=-1

        #pass by pass data
        self.t_hist=np.zeros(beamtim.N_pass_tot+1,float)
        self.nel_impact_hist_tot=np.zeros((beamtim.N_pass_tot+1,impact_man.Nxg_hist),float)
        self.nel_impact_hist_scrub=np.zeros((beamtim.N_pass_tot+1,impact_man.Nxg_hist),float)
        self.energ_eV_impact_hist=np.zeros((beamtim.N_pass_tot+1,impact_man.Nxg_hist),float)
        self.nel_hist_line=np.zeros(impact_man.Nxg_hist,float)
        self.nel_hist=np.zeros((beamtim.N_pass_tot+1,impact_man.Nxg_hist),float)
        self.N_mp_impact_pass=np.zeros(beamtim.N_pass_tot+1)
        self.N_mp_corrected_pass=np.zeros(beamtim.N_pass_tot+1)
        self.N_mp_pass=np.zeros(beamtim.N_pass_tot+1)
        self.N_mp_ref_pass=np.zeros(beamtim.N_pass_tot+1)

        if impact_man.flag_seg:
                self.nel_hist_impact_seg=np.zeros((beamtim.N_pass_tot+1,impact_man.chamb.N_vert),float)
                self.energ_eV_impact_seg=np.zeros((beamtim.N_pass_tot+1,impact_man.chamb.N_vert),float)
        else:
                self.nel_hist_impact_seg=-1
                self.energ_eV_impact_seg=-1


        #logfile and progress file
        self.logfile_path=logfile_path
        self.progress_path=progress_path


        if flag_presence_sec_beams:
            self.t_sec_beams = beamtim.t[::dec_fac_secbeam_prof]
            self.sec_beam_profiles = []
            for sec_beam in sec_beams_list:
                self.sec_beam_profiles.append(sec_beam.lam_t_array[::dec_fac_secbeam_prof])
            self.sec_beam_profiles=np.array(self.sec_beam_profiles)
        else:
            self.t_sec_beams = -1
            self.sec_beam_profiles = -1

        # initialize electron density probes
        self.flag_el_dens_probes = False
        self.x_el_dens_probes = -1
        self.y_el_dens_probes = -1
        self.r_el_dens_probes = -1
        self.el_dens_at_probes = -1
        if len(el_density_probes)>0:
            self.flag_el_dens_probes = True
            self.N_el_dens_probes = len(el_density_probes)
            self.el_dens_at_probes = np.zeros((self.N_el_dens_probes, len(self.t_dec)))
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



        self.stopfile = stopfile

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
            self.nel_hist_det=np.zeros((beamtim.N_pass_tot+1,self.Nxg_hist_det),float)



        print('Done init pyecloud_saver.')
        flog=open(self.logfile_path,'a')
        timestr = time.strftime("%d %b %Y %H:%M:%S", time.localtime())
        flog.write('Initialization finished on %s\n'%timestr)
        flog.close()

    def witness(self, MP_e, beamtim, spacech_ele, impact_man,
        dynamics,gas_ion_flag,resgasion,t_ion,
        t_sc_ON, photoem_flag, phemiss,flag_presence_sec_beams,sec_beams_list):

        #MP state save
        if self.flag_save_MP_state:
            if  (MP_e.N_mp>0) and (self.i_obs<self.N_obs):
                if (beamtim.tt_curr>=self.t_obs[self.i_obs]):
                    filename_MP_state='MP_state_%d'%(self.i_obs)
                    sio.savemat(filename_MP_state,{'tt':beamtim.tt_curr,'N_mp':MP_e.N_mp, 'x_mp':MP_e.x_mp[0:MP_e.N_mp], 'y_mp':MP_e.y_mp[0:MP_e.N_mp], 'z_mp':MP_e.z_mp[0:MP_e.N_mp],\
                                                    'vx_mp':MP_e.vx_mp[0:MP_e.N_mp], 'vy_mp':MP_e.vy_mp[0:MP_e.N_mp], 'vz_mp':MP_e.vz_mp[0:MP_e.N_mp], 'nel_mp':MP_e.nel_mp[0:MP_e.N_mp]},oned_as='row')

                    print('Save MP state on ' + filename_MP_state)
                    self.i_obs=self.i_obs+1

        #Simulation state save
        if self.flag_save_simulation_state:
            if  self.i_obs_sim<self.N_obs_sim:
                if (beamtim.tt_curr>=self.t_obs_sim[self.i_obs_sim]):
                    filename_simulation_state='simulation_state_%d.pkl'%(self.i_obs_sim)

                    temp_luobj = spacech_ele.PyPICobj.luobj
                    spacech_ele.luobj=None
                    spacech_ele.PyPICobj.luobj=None

                    #~ dynamics.get_B=None

                    dict_state = {
                    'beamtim':beamtim,
                    'MP_e':MP_e,
                    'dynamics':dynamics,
                    'impact_man':impact_man,
                    'gas_ion_flag':gas_ion_flag,
                    'resgasion':resgasion,
                    't_ion':t_ion,
                    'spacech_ele':spacech_ele,
                    't_sc_ON':t_sc_ON,
                    'photoem_flag':photoem_flag,
                    'phemiss':phemiss,
                    'flag_presence_sec_beams':flag_presence_sec_beams,
                    'sec_beams_list':sec_beams_list}



                    with open(filename_simulation_state, 'wb') as fid:
                        # use best protocol available
                        pickle.dump(dict_state, fid, protocol=-1)

                    spacech_ele.PyPICobj.luobj = temp_luobj

                    print('Save simulation state on ' + filename_simulation_state)
                    self.i_obs_sim=self.i_obs_sim+1

        # Energy histogram saver
        if (np.mod(beamtim.ii_curr,self.Nst_En_hist)==0):
            self.En_hist[self.i_En_hist,:]=impact_man.En_hist_line
            self.t_En_hist[self.i_En_hist]=beamtim.tt_curr
            impact_man.reset_En_hist_line()
            self.i_En_hist=self.i_En_hist+1

            if self.flag_cos_angle_hist:
                self.cos_angle_hist[self.i_En_hist,:] = impact_man.cos_angle_hist
                impact_man.reset_cos_angle_hist()

        #Space charge electrostatic energy
        if spacech_ele.flag_recomputed_sc:
            self.t_sc_video.append(beamtim.tt_curr)
            self.U_sc_eV.append(spacech_ele.U_sc_eV_stp)

        #save rho video
        if self.flag_video:
            if not os.path.exists('rho_video'):
                os.makedirs('rho_video')
            if self.rho_video is None:
                self.rho_video=[]
                self.t_video=[]
            if spacech_ele.flag_recomputed_sc:
                self.rho_video.append(spacech_ele.rho)
                self.t_video.append(beamtim.tt_curr)
            if beamtim.flag_new_bunch_pass:
                self.rho_video=np.array(self.rho_video)
                self.t_video=np.array(self.t_video)
                filename_rho='rho_video/rho_pass%d.mat'%(beamtim.pass_numb-1)
                print('Saving %s'%filename_rho)
                sio.savemat(filename_rho,{'xg_sc':spacech_ele.xg,'yg_sc':spacech_ele.yg,'t_video':self.t_video,'rho_video':self.rho_video},oned_as='row')
                print('Done')
                self.rho_video=[]
                self.t_video=[]

        #save efield video
        if self.flag_sc_video:
            if not os.path.exists('efield_video'):
                os.makedirs('efield_video')
            if self.efx_video is None:
                self.efx_video=[]
                self.efy_video=[]
                self.t_efield_video=[]
            if spacech_ele.flag_recomputed_sc:
                self.efx_video.append(spacech_ele.efx)
                self.efy_video.append(spacech_ele.efy)
                self.t_efield_video.append(beamtim.tt_curr)
            if beamtim.flag_new_bunch_pass:
                self.efx_video=np.array(self.efx_video)
                self.efy_video=np.array(self.efy_video)
                self.t_efield_video=np.array(self.t_efield_video)
                filename_efield='efield_video/efield_pass%d.mat'%(beamtim.pass_numb-1)
                print('Saving %s'%filename_efield)
                sio.savemat(filename_efield,{'xg_sc':spacech_ele.xg,'yg_sc':spacech_ele.yg,'t_efield_video':self.t_efield_video,
                                          'efx_video':self.efx_video, 'efy_video':self.efy_video},oned_as='row')
                print('Done')
                self.efx_video=[]
                self.efy_video=[]
                self.t_efield_video=[]


        #save step by step data
        # Vars to be accumulated
        self.Nel_impact_last_step_group += impact_man.Nel_impact_last_step
        self.Nel_emit_last_step_group += impact_man.Nel_emit_last_step
        self.En_imp_last_step_group_eV += impact_man.En_imp_last_step_eV
        self.En_emit_last_step_group_eV += impact_man.En_emit_last_step_eV


        if np.mod(beamtim.ii_curr, self.dec_fact_out)==0:
            ii_curr_dec=beamtim.ii_curr/self.dec_fact_out
            self.Nel_imp_time[ii_curr_dec] = self.Nel_impact_last_step_group
            self.Nel_emit_time[ii_curr_dec] = self.Nel_emit_last_step_group
            self.En_imp_eV_time[ii_curr_dec] = self.En_imp_last_step_group_eV
            self.En_emit_eV_time[ii_curr_dec] = self.En_emit_last_step_group_eV

            self.Nel_impact_last_step_group = 0
            self.Nel_emit_last_step_group = 0
            self.En_imp_last_step_group_eV = 0
            self.En_emit_last_step_group_eV = 0


            self.Nel_time[ii_curr_dec]=np.sum(MP_e.nel_mp[0:MP_e.N_mp]);
            self.En_kin_eV_time[ii_curr_dec]=np.sum(0.5/qm*MP_e.nel_mp[0:MP_e.N_mp]*(MP_e.vx_mp[0:MP_e.N_mp]*MP_e.vx_mp[0:MP_e.N_mp]+MP_e.vy_mp[0:MP_e.N_mp]*MP_e.vy_mp[0:MP_e.N_mp]+MP_e.vz_mp[0:MP_e.N_mp]*MP_e.vz_mp[0:MP_e.N_mp]));

            flag_center=((MP_e.x_mp**2 + MP_e.y_mp**2)<self.r_center**2);
            flag_center[MP_e.N_mp:]=False
            self.cen_density[ii_curr_dec]=np.sum(MP_e.nel_mp[flag_center])/(np.pi*self.r_center*self.r_center)

            if self.flag_el_dens_probes:
                for ii in xrange(self.N_el_dens_probes):
                    flag_center=((MP_e.x_mp-self.x_el_dens_probes[ii])**2 + (MP_e.y_mp-self.y_el_dens_probes[ii])**2)<self.r_el_dens_probes[ii]**2;
                    flag_center[MP_e.N_mp:]=False
                    self.el_dens_at_probes[ii, ii_curr_dec]=np.sum(MP_e.nel_mp[flag_center])/(np.pi*self.r_el_dens_probes[ii]**2)


            if self.flag_detailed_MP_info==1:
                self.N_mp_time[ii_curr_dec]=MP_e.N_mp


        if beamtim.flag_new_bunch_pass:

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


            self.nel_hist[beamtim.pass_numb,:]=self.nel_hist_line
            self.t_hist[beamtim.pass_numb]=beamtim.tt_curr
            self.nel_impact_hist_tot[beamtim.pass_numb,:]=impact_man.nel_impact_hist_tot
            impact_man.reset_impact_hist_tot()
            self.nel_impact_hist_scrub[beamtim.pass_numb,:]=impact_man.nel_impact_hist_scrub
            impact_man.reset_impact_hist_scrub()
            self.energ_eV_impact_hist[beamtim.pass_numb,:]=impact_man.energ_eV_impact_hist
            impact_man.reset_energ_eV_impact_hist()

            self.N_mp_impact_pass[beamtim.pass_numb]=impact_man.chamb.N_mp_impact
            self.N_mp_corrected_pass[beamtim.pass_numb]=impact_man.chamb.N_mp_corrected
            self.N_mp_pass[beamtim.pass_numb]=MP_e.N_mp
            self.N_mp_ref_pass[beamtim.pass_numb]=MP_e.nel_mp_ref

            if impact_man.flag_seg:
                self.nel_hist_impact_seg[beamtim.pass_numb,:]=impact_man.nel_hist_impact_seg
                impact_man.reset_hist_impact_seg()

            if impact_man.flag_seg:
                self.energ_eV_impact_seg[beamtim.pass_numb,:]=impact_man.energ_eV_impact_seg
                impact_man.reset_energ_impact_seg()

            if self.flag_hist_det:
                self.nel_hist_det[beamtim.pass_numb,:]=self.nel_hist_det_line

            #self.t_dec = beamtim.t[::dec_fact_out]
            #self.lam_t_array_dec = beamtim.lam_t_array[::dec_fact_out]

            saved_dict = {
                    'Nel_timep':            self.Nel_time,
                    't':                    self.t_dec,
                    't_hist':               self.t_hist,
                    'nel_hist':             self.nel_hist,
                    'xg_hist':              impact_man.xg_hist,
                    'nel_impact_hist_tot':  self.nel_impact_hist_tot,
                    'nel_impact_hist_scrub':self.nel_impact_hist_scrub,
                    'energ_eV_impact_hist': self.energ_eV_impact_hist,
                    'cen_density':          self.cen_density,
                    'Nel_imp_time':         self.Nel_imp_time,
                    'Nel_emit_time':        self.Nel_emit_time,
                    'En_imp_eV_time':       self.En_imp_eV_time,
                    'En_emit_eV_time':      self.En_emit_eV_time,
                    'En_g_hist':            impact_man.En_g_hist,
                    'En_hist':              self.En_hist,
                    't_En_hist':            self.t_En_hist,
                    'lam_t_array':          self.lam_t_array_dec,
                    'b_spac':               beamtim.b_spac,
                    't_sc_video':           np.array(self.t_sc_video),
                    'U_sc_eV':              np.array(self.U_sc_eV),
                    'En_kin_eV_time':       self.En_kin_eV_time,
                    'N_mp_impact_pass':     self.N_mp_impact_pass,
                    'N_mp_corrected_pass':  self.N_mp_corrected_pass,
                    'N_mp_pass':            self.N_mp_pass,
                    'N_mp_time':            self.N_mp_time,
                    'N_mp_ref_pass':        self.N_mp_ref_pass,
                    'nel_hist_impact_seg':  self.nel_hist_impact_seg,
                    'energ_eV_impact_seg':  self.energ_eV_impact_seg,
                    't_sec_beams':          self.t_sec_beams,
                    'sec_beam_profiles':    self.sec_beam_profiles,
                    'el_dens_at_probes':    self.el_dens_at_probes,
                    'x_el_dens_probes':     self.x_el_dens_probes,
                    'y_el_dens_probes':     self.y_el_dens_probes,
                    'r_el_dens_probes':     self.r_el_dens_probes,
                    'nel_hist_det':         self.nel_hist_det,
                    'xg_hist_det':          self.xg_hist_det,
                    'dec_fact_out':         self.dec_fact_out,
                }
            if self.flag_cos_angle_hist:
                    saved_dict['cos_angle_hist'] = self.cos_angle_hist
                    saved_dict['xg_hist_cos_angle'] = self.xg_hist_cos_angle

            sio.savemat(self.filen_main_outp, saved_dict, oned_as='row')


            # logfile and progressfile
            timestr = time.strftime("%d %b %Y %H:%M:%S", time.localtime())

            string_tolog= timestr+(' pass. %d/%d Nel_tot=%e N_mp=%d\n'%(beamtim.pass_numb,beamtim.N_pass_tot,np.sum(MP_e.nel_mp[0:MP_e.N_mp]),MP_e.N_mp))


#             flog=open(self.logfile_path,'a')
#             flog.write(string_tolog)
#             flog.close()
            try:
                with open(self.logfile_path,'a') as flog:
                    flog.write(string_tolog)
            except IOError as err:
                print('Got: ',err)
                print('while trying to write the following line on logfile:')
                print(string_tolog)

#             flog=open(self.progress_path,'w')
#             flog.write(('%f'%(float(beamtim.pass_numb)/float(beamtim.N_pass_tot))))
#             flog.close()
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

        return impact_man

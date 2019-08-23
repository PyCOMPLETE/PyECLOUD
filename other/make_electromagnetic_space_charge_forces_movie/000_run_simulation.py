import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys, os
BIN = os.path.expanduser("../../../")
sys.path.append(BIN)
sim_input_folder = '../../testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_em_tracking'
from PyECLOUD.buildup_simulation import BuildupSimulation
from scipy.constants import e as qe
import PyECLOUD.hist_for as histf
import PyECLOUD.myloadmat_to_obj as mlo
import scipy.io as sio
import PyECLOUD.rhocompute as rhocom

def _forces_movie_init(sim,sim_input_folder):
    sim.cloud_list[0].pyeclsaver.Fefx_video = None
    sim.cloud_list[0].pyeclsaver.Fefy_video = None
    sim.cloud_list[0].pyeclsaver.Fbfx_video = None
    sim.cloud_list[0].pyeclsaver.Fbfy_video = None
    sim.cloud_list[0].pyeclsaver.Fbfz_video = None
    sim.cloud_list[0].pyeclsaver.t_efield_video_new = None
    sim.cloud_list[0].pyeclsaver.sim_input_folder = sim_input_folder
    sim.cloud_list[0].pyeclsaver._rho_video_init(True)
    return 0


def _forces_movie_save(sim):
    spacech_ele = sim.spacech_ele
    MP_e = sim.cloud_list[0].MP_e
    sim_input_folder = sim.cloud_list[0].pyeclsaver.sim_input_folder
    beamtim = sim.beamtim
    pyeclsaver = sim.cloud_list[0].pyeclsaver
    #get the fields (on the particles)
    efx, efy, bfx, bfy, bfz = spacech_ele.get_sc_em_field(MP_e)
    #compute the forces (on the macro-particles)
    Fefx_mp = qe*np.multiply(MP_e.nel_mp[0:MP_e.N_mp],efx)
    Fefy_mp = qe*np.multiply(MP_e.nel_mp[0:MP_e.N_mp],efy)
    Fbfx_mp = qe*np.multiply(MP_e.nel_mp[0:MP_e.N_mp],np.multiply(MP_e.vy_mp[0:MP_e.N_mp],bfz)-qe*np.multiply(MP_e.vz_mp[0:MP_e.N_mp],bfy))
    Fbfy_mp = qe*np.multiply(MP_e.nel_mp[0:MP_e.N_mp],np.multiply(MP_e.vz_mp[0:MP_e.N_mp],bfx)-qe*np.multiply(MP_e.vx_mp[0:MP_e.N_mp],bfz))
    Fbfz_mp = qe*np.multiply(MP_e.nel_mp[0:MP_e.N_mp],np.multiply(MP_e.vx_mp[0:MP_e.N_mp],bfy)-qe*np.multiply(MP_e.vy_mp[0:MP_e.N_mp],bfx))
    #interpolate the forces to the mesh nodes
    Fefx = rhocom.compute_sc_rho(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp],Fefx_mp,spacech_ele.bias_x,spacech_ele.bias_y,spacech_ele.Dh,spacech_ele.Nxg,spacech_ele.Nyg)
    Fefy = rhocom.compute_sc_rho(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp],Fefy_mp,spacech_ele.bias_x,spacech_ele.bias_y,spacech_ele.Dh,spacech_ele.Nxg,spacech_ele.Nyg)
    Fbfx = rhocom.compute_sc_rho(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp],Fbfx_mp,spacech_ele.bias_x,spacech_ele.bias_y,spacech_ele.Dh,spacech_ele.Nxg,spacech_ele.Nyg)
    Fbfy = rhocom.compute_sc_rho(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp],Fbfy_mp,spacech_ele.bias_x,spacech_ele.bias_y,spacech_ele.Dh,spacech_ele.Nxg,spacech_ele.Nyg)
    Fbfz = rhocom.compute_sc_rho(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp],Fbfz_mp,spacech_ele.bias_x,spacech_ele.bias_y,spacech_ele.Dh,spacech_ele.Nxg,spacech_ele.Nyg)
    #save sc forces video
    if not os.path.exists(sim_input_folder+ '/sc_forces'):
        os.makedirs(sim_input_folder+ '/sc_forces')
    if pyeclsaver.Fefx_video is None:
        pyeclsaver.Fefx_video = []
        pyeclsaver.Fefy_video = []
        pyeclsaver.Fbfx_video = []
        pyeclsaver.Fbfy_video = []
        pyeclsaver.Fbfz_video = []
        pyeclsaver.t_efield_video_new = []
    if spacech_ele.last_recomputation_check:
        pyeclsaver.Fefx_video.append(Fefx)
        pyeclsaver.Fefy_video.append(Fefy)
        pyeclsaver.Fbfx_video.append(Fbfx)
        pyeclsaver.Fbfy_video.append(Fbfy)
        pyeclsaver.Fbfz_video.append(Fbfz)
        pyeclsaver.t_efield_video_new.append(beamtim.tt_curr)
    if beamtim.flag_new_bunch_pass:
        pyeclsaver.Fefx_video = np.array(pyeclsaver.Fefx_video)
        pyeclsaver.Fefy_video = np.array(pyeclsaver.Fefy_video)
        pyeclsaver.Fbfx_video = np.array(pyeclsaver.Fbfx_video)
        pyeclsaver.Fbfy_video = np.array(pyeclsaver.Fbfy_video)
        pyeclsaver.Fbfz_video = np.array(pyeclsaver.Fbfz_video)
        pyeclsaver.t_efield_video_new = np.array(pyeclsaver.t_efield_video_new)
        filename_sc_forces = pyeclsaver.folder_outp + '/sc_forces/sc_forces_pass%d.mat'%(beamtim.pass_numb - 1)
        print('Saving %s'%filename_sc_forces)
        sio.savemat(filename_sc_forces, {'xg_sc': spacech_ele.xg, 'yg_sc': spacech_ele.yg, 't_efield_video_new': pyeclsaver.t_efield_video_new,
                                  'Fefx_video': pyeclsaver.Fefx_video, 'Fefy_video': pyeclsaver.Fefy_video, 'Fbfx_video': pyeclsaver.Fbfx_video,
                                  'Fbfy_video': pyeclsaver.Fbfy_video, 'Fbfz_video': pyeclsaver.Fbfz_video}, oned_as='row')
        print('Done')
        pyeclsaver.Fefx_video = []
        pyeclsaver.Fefy_video = []
        pyeclsaver.Fbfx_video = []
        pyeclsaver.Fbfy_video = []
        pyeclsaver.Fbfz_video = []
        pyeclsaver.t_efield_video_new = []
    return 0


step_by_step_custom_observables = {
        'dummy2': lambda sim : _forces_movie_save(sim)
        }


sim = BuildupSimulation(pyecl_input_folder=sim_input_folder,
        filen_main_outp='./Pyecltest.mat',
        extract_sey=False,
        step_by_step_custom_observables = step_by_step_custom_observables)

_forces_movie_init(sim,sim_input_folder)

sim.run()

ob1 = mlo.myloadmat('Pyecltest.mat')
cloud = sim.cloud_list[0]
pyeclsaver = cloud.pyeclsaver
ob1['forces_ratio_hist'] = pyeclsaver.forces_ratio_hist
ob1['forces_ratio_hist_bins'] = pyeclsaver.forces_ratio_hist_bins
ob1['t_forces_ratio_hist'] = pyeclsaver.t_forces_ratio_hist
sio.savemat('Pyecltest_enhanced.mat', ob1, oned_as='row')

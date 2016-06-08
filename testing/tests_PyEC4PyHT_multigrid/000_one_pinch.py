import sys, os
BIN=os.path.expanduser('../../../')
sys.path.insert(0,BIN)

from scipy.constants import c, e, m_p
import numpy as np
import pylab as pl				
import myfilemanager as mlm
import PyECLOUD.mystyle as ms

n_turns = 512
n_segments =1
machine_configuration = '6.5_TeV_collision_tunes'

z_cut = 2.5e-9*c
n_slices = 150
L_ecloud = 1000.

sigma_z = 10e-2
epsn_x = 2.5e-6
epsn_y = 2.5e-6



# Define the machine
#============================
from LHC import LHC
machine = LHC(machine_configuration=machine_configuration,
                        optics_mode='smooth', n_segments=n_segments)

bunch = machine.generate_6D_Gaussian_bunch(
                                        n_macroparticles=3000000, intensity=1e11,
                                        epsn_x=epsn_x, epsn_y=epsn_y, sigma_z=sigma_z)


             
from PyHEADTAIL.particles.slicing import UniformBinSlicer
slicer = UniformBinSlicer(n_slices = n_slices, z_cuts=(-z_cut, z_cut) )

x_beam_offset = 0.
y_beam_offset = 0.
D_probe = bunch.sigma_x()/2

probes_position = [{'x' : x_beam_offset, 'y': y_beam_offset+D_probe},
                    {'x' : x_beam_offset, 'y': y_beam_offset-D_probe},
                    {'x' : x_beam_offset+D_probe, 'y': y_beam_offset},
                    {'x' : x_beam_offset-D_probe, 'y': y_beam_offset},
                    {'x' : x_beam_offset, 'y': y_beam_offset+(2*D_probe)},
                    {'x' : x_beam_offset, 'y': y_beam_offset-(2*D_probe)},
                    {'x' : x_beam_offset+(2*D_probe), 'y': y_beam_offset},
                    {'x' : x_beam_offset-(2*D_probe), 'y': y_beam_offset}]

                
import PyECLOUD.PyEC4PyHT as PyEC4PyHT                        
ecloud_singlegrid = PyEC4PyHT.Ecloud(
        L_ecloud=L_ecloud, slicer=slicer,
        Dt_ref=20e-12, pyecl_input_folder='./pyecloud_config',
        chamb_type = 'polyg' ,
        filename_chm= 'LHC_chm_ver.mat', Dh_sc = 0.1*bunch.sigma_x(),
        init_unif_edens_flag=1,
        init_unif_edens=1e7,
        N_mp_max = 3000000,
        nel_mp_ref_0 = 1e7/(0.7*3000000),
        B_multip = [0.],
        x_beam_offset = x_beam_offset,
        y_beam_offset = y_beam_offset,
        probes_position = probes_position)
        
ecloud_multigrid = PyEC4PyHT.Ecloud(
        L_ecloud=L_ecloud, slicer=slicer,
        Dt_ref=20e-12, pyecl_input_folder='./pyecloud_config',
        chamb_type = 'polyg' ,
        filename_chm= 'LHC_chm_ver.mat', Dh_sc=1e-3,
        init_unif_edens_flag=1,
        init_unif_edens=1e7,
        N_mp_max = 3000000,
        nel_mp_ref_0 = 1e7/(0.7*3000000),
        B_multip = [0.],
        PyPICmode = 'ShortleyWeller_WithTelescopicGrids',
        f_telescope = 0.5,
        target_grid = {'x_min_target':-5*bunch.sigma_x(), 'x_max_target':5*bunch.sigma_x(),'y_min_target':-5*bunch.sigma_y(),'y_max_target':5*bunch.sigma_y(),'Dh_target':0.1*bunch.sigma_x()},
        N_nodes_discard = 3.,
        N_min_Dh_main = 10,
        x_beam_offset = x_beam_offset,
        y_beam_offset = y_beam_offset,
        probes_position = probes_position)
        

                                        								
#~ new_ecloud.save_ele_distributions_last_track = True
#~ new_ecloud.save_ele_potential_and_field = True
#~ new_ecloud.track_once_and_replace_with_recorded_field_map(bunch_for_field_calculation, delete_ecloud_data=False)
#~ machine.one_turn_map.append(new_ecloud)

import time
N_rep_test = 2
print 'Ecloud track %d times'%N_rep_test
t_start_sw = time.mktime(time.localtime())
for _ in xrange(N_rep_test):
    ecloud_singlegrid.track(bunch) 
t_stop_sw = time.mktime(time.localtime())
t_sw_single = t_stop_sw-t_start_sw
print 'Singlegrid tracking time ', t_sw_single /N_rep_test 

import time
N_rep_test = 2
print 'Ecloud track %d times'%N_rep_test
t_start_sw = time.mktime(time.localtime())
for _ in xrange(N_rep_test):
    ecloud_multigrid.track(bunch) 
t_stop_sw = time.mktime(time.localtime())
t_sw_single = t_stop_sw-t_start_sw
print 'Multigrid tracking time ', t_sw_single /N_rep_test 


                     
flag_plot_efield_at_probes = True
if flag_plot_efield_at_probes:

        
        slices = bunch.get_slices(ecloud_singlegrid.slicer)
            
        pl.close('all')
        
        import matplotlib.gridspec as gridspec
        fig = pl.figure(1, figsize=(7,8))
        fig.patch.set_facecolor('w')
        gs1 = gridspec.GridSpec(1, 1)
        gs2 = gridspec.GridSpec(2, 1)
        
        sp1 = fig.add_subplot(gs1[0])
        obcham =  mlm.myloadmat_to_obj( 'LHC_chm_ver.mat') 	
        sp1.plot(obcham.Vx*1e3, obcham.Vy*1e3, 'b')
        sp1.plot(ecloud_singlegrid.x_probes*1e3, ecloud_singlegrid.y_probes*1e3, 'or', markersize=3)
        sp1.plot(x_beam_offset, y_beam_offset, '*k', markersize=4)
        #~ sp1.set_xlim([-0.020, 0.020])
        #~ sp1.set_ylim([-0.020, 0.020])
        sp1.set_ylabel('y [mm]')
        sp1.set_xlabel('x [mm]')
        sp1.axis('equal')
        ms.sciy()
        ms.scix()
        #sp1.axis('equal')
        sp1.grid('on')
        
        
        sp2 = fig.add_subplot(gs2[0])
        sp2.plot(slices.z_centers*1e2, ecloud_singlegrid.Ex_ele_last_track_at_probes, 'b.', markersize=6)
        sp2.plot(slices.z_centers*1e2, ecloud_multigrid.Ex_ele_last_track_at_probes, 'r.', markersize=6)
        sp2.set_ylabel('Ex at probes [V/m]')
        sp2.set_xlabel('z [cm]')
        sp2.grid('on')
        ms.sciy()
        
        sp3 = fig.add_subplot(gs2[1])
        sp3.plot(slices.z_centers*1e2, ecloud_singlegrid.Ey_ele_last_track_at_probes, 'b.', markersize=6)
        sp3.plot(slices.z_centers*1e2, ecloud_multigrid.Ey_ele_last_track_at_probes, 'r.', markersize=6)

        sp3.set_ylabel('Ey at probes [V/m]')
        sp3.set_xlabel('z [cm]')
        sp3.grid('on')
        
        ms.sciy()

        gs1.tight_layout(fig, rect=[0.25, 0.6, 0.75, 0.98], pad=1.08, h_pad=.9)     #[left, bottom, right, top]
        gs2.tight_layout(fig, rect=[0, 0, 1, 0.65], pad=1.08, h_pad=.9)     #[left, bottom, right, top]
            
        pl.savefig('1.png', dpi=300)
        pl.show()
        
        #~ dict_to_mat = {\
        #~ 'n_slices':slices.n_slices,\
        #~ 'z_centers':slices.z_centers,\
        #~ 'x_probes':new_ecloud.x_probes,\
        #~ 'y_probes':new_ecloud.y_probes,\
        #~ 'Ex_ele_last_track_at_probes':new_ecloud.Ex_ele_last_track_at_probes,\
        #~ 'Ey_ele_last_track_at_probes':new_ecloud.Ey_ele_last_track_at_probes,\
        #~ 'x_beam_offset':x_beam_offset,\
        #~ 'y_beam_offset':y_beam_offset}

        #~ import scipy.io as sio
        #~ sio.savemat('var_for_analysis_ref', dict_to_mat, oned_as='row')


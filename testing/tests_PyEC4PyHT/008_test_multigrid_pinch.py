import sys, os
BIN=os.path.expanduser('../../../')
sys.path.insert(0,BIN)

from scipy.constants import c, e, m_p
import numpy as np
import pylab as pl				
import myfilemanager as mlm
import PyECLOUD.mystyle as ms

n_segments =1
machine_configuration = '6.5_TeV_collision_tunes'
p0_GeV = 450

z_cut = 2.5e-9*c
n_slices = 150
L_ecloud = 1000.

sigma_z = 10e-2
epsn_x = 2.5e-6
epsn_y = 2.5e-6

sparse_solver = 'PyKLU'#'scipy_slu'


# Define the machine
#============================
from machines_for_testing import LHC
machine = LHC(machine_configuration=machine_configuration,
                        optics_mode='smooth', n_segments=n_segments, p0=p0_GeV*1e9*e/c)

bunch = machine.generate_6D_Gaussian_bunch(
                                        n_macroparticles=3000000, intensity=1e11,
                                        epsn_x=epsn_x, epsn_y=epsn_y, sigma_z=sigma_z)

             
from PyHEADTAIL.particles.slicing import UniformBinSlicer
slicer = UniformBinSlicer(n_slices = n_slices, z_cuts=(-z_cut, z_cut) )

x_beam_offset = 0.
y_beam_offset = 0.
Dx_probe = bunch.sigma_x()
Dy_probe = bunch.sigma_y()

probes_position = [{'x' : x_beam_offset, 'y': y_beam_offset+Dy_probe},
                    {'x' : x_beam_offset, 'y': y_beam_offset-Dy_probe},
                    {'x' : x_beam_offset+Dx_probe, 'y': y_beam_offset},
                    {'x' : x_beam_offset-Dx_probe, 'y': y_beam_offset},
                    {'x' : x_beam_offset, 'y': y_beam_offset+(2*Dy_probe)},
                    {'x' : x_beam_offset, 'y': y_beam_offset-(2*Dy_probe)},
                    {'x' : x_beam_offset+(2*Dx_probe), 'y': y_beam_offset},
                    {'x' : x_beam_offset-(2*Dx_probe), 'y': y_beam_offset}]
                    
                
import PyECLOUD.PyEC4PyHT as PyEC4PyHT                        
ecloud_singlegrid = PyEC4PyHT.Ecloud(
        L_ecloud=L_ecloud, slicer=slicer,
        Dt_ref=20e-12, pyecl_input_folder='./pyecloud_config_LHC',
        chamb_type = 'polyg' ,
        filename_chm= 'LHC_chm_ver.mat', Dh_sc = .2*bunch.sigma_x(),
        init_unif_edens_flag=1,
        init_unif_edens=1e7,
        N_mp_max = 3000000,
        nel_mp_ref_0 = 1e7/(0.7*3000000),
        B_multip = [0.],
        x_beam_offset = x_beam_offset,
        y_beam_offset = y_beam_offset,
        probes_position = probes_position,
        sparse_solver = sparse_solver)
        
ecloud_multigrid = PyEC4PyHT.Ecloud(
        L_ecloud=L_ecloud, slicer=slicer,
        Dt_ref=20e-12, pyecl_input_folder='./pyecloud_config_LHC',
        chamb_type = 'polyg' ,
        filename_chm= 'LHC_chm_ver.mat', Dh_sc=1e-3,
        init_unif_edens_flag=1,
        init_unif_edens=1e7,
        N_mp_max = 3000000,
        nel_mp_ref_0 = 1e7/(0.7*3000000),
        B_multip = [0.],
        PyPICmode = 'ShortleyWeller_WithTelescopicGrids',
        f_telescope = 0.3,
        target_grid = {'x_min_target':-5*bunch.sigma_x(), 'x_max_target':5*bunch.sigma_x(),
                       'y_min_target':-5*bunch.sigma_y(),'y_max_target':5*bunch.sigma_y(),
                       'Dh_target':.2*bunch.sigma_x()},
        N_nodes_discard = 10.,
        N_min_Dh_main = 10,
        x_beam_offset = x_beam_offset,
        y_beam_offset = y_beam_offset,
        probes_position = probes_position,
        sparse_solver = sparse_solver)
     
     
# Singlegrid tracking time 
import time
N_rep_test = 1
print 'Ecloud track %d times'%N_rep_test
t_start_sw = time.mktime(time.localtime())
for _ in xrange(N_rep_test):
    ecloud_singlegrid.track(bunch) 
t_stop_sw = time.mktime(time.localtime())
t_sw_single = t_stop_sw-t_start_sw
print 'Singlegrid tracking time ', t_sw_single /N_rep_test 
# Multigrid tracking time 
print 'Ecloud track %d times'%N_rep_test
t_start_sw = time.mktime(time.localtime())
for _ in xrange(N_rep_test):
    ecloud_multigrid.track(bunch) 
t_stop_sw = time.mktime(time.localtime())
t_sw_single = t_stop_sw-t_start_sw
print 'Multigrid tracking time ', t_sw_single /N_rep_test 


                    
slices = bunch.get_slices(ecloud_singlegrid.slicer)


pl.close('all')

pic_singlegrid = ecloud_singlegrid.spacech_ele.PyPICobj
pic_multigrid = ecloud_multigrid.spacech_ele.PyPICobj

# Plot field at probes 
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
for ii in xrange(pic_multigrid.n_grids):
    sp1.plot(pic_multigrid.pic_list[ii].pic_internal.chamb.Vx*1e3, pic_multigrid.pic_list[ii].pic_internal.chamb.Vy*1e3, '.-k')
sp1.set_ylabel('y [mm]')
sp1.set_xlabel('x [mm]')
sp1.axis('equal')
ms.sciy()
ms.scix()
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



########################
#plot field maps last slice
########################    

#probes for map
Dh_test = .5e-4
x_grid_probes = np.arange(np.min(pic_singlegrid.xg), np.max(pic_singlegrid.xg)+Dh_test, Dh_test)
y_grid_probes = np.arange(np.min(pic_singlegrid.yg), np.max(pic_singlegrid.yg), Dh_test)

[xn, yn]=np.meshgrid(x_grid_probes,y_grid_probes)
xn=xn.T
xn=xn.flatten()
yn=yn.T
yn=yn.flatten()

#pic gather field 
Ex_singlegrid_n, Ey_singlegrid_n = pic_singlegrid.gather(xn, yn)	
Ex_singlegrid_matrix=np.reshape(Ex_singlegrid_n,(len(y_grid_probes),len(x_grid_probes)), 'F').T
Ey_singlegrid_matrix=np.reshape(Ey_singlegrid_n,(len(y_grid_probes),len(x_grid_probes)), 'F').T

Ex_multigrid_n, Ey_multigrid_n = pic_multigrid.gather(xn, yn)	
Ex_multigrid_matrix=np.reshape(Ex_multigrid_n,(len(y_grid_probes),len(x_grid_probes)), 'F').T
Ey_multigrid_matrix=np.reshape(Ey_multigrid_n,(len(y_grid_probes),len(x_grid_probes)), 'F').T

#pic gather rho 
rho_multigrid_n = pic_multigrid.gather_rho(xn, yn)	
rho_multigrid_matrix=np.reshape(rho_multigrid_n,(len(y_grid_probes),len(x_grid_probes)), 'F').T

rho_singlegrid_n = pic_singlegrid.gather_rho(xn, yn)	
rho_singlegrid_matrix=np.reshape(rho_singlegrid_n,(len(y_grid_probes),len(x_grid_probes)), 'F').T

#plot
vmin = -7; vmax = -2
pl.figure(2, figsize=(12, 6)).patch.set_facecolor('w')
pl.subplots_adjust(hspace = 0.5, wspace = 0.3)
sp1 = pl.subplot(2,2,1)
pl.pcolormesh(x_grid_probes, y_grid_probes, 
	np.log10(np.sqrt(Ex_singlegrid_matrix**2+Ey_singlegrid_matrix**2).T), vmin=vmin, vmax=vmax)
for ii in xrange(pic_multigrid.n_grids):
    sp1.plot(pic_multigrid.pic_list[ii].pic_internal.chamb.Vx, pic_multigrid.pic_list[ii].pic_internal.chamb.Vy, '.-')
pl.xlabel('x [m]')
pl.ylabel('y [m]')
cb=pl.colorbar(); pl.axis('equal')
cb.formatter.set_powerlimits((0, 0))
cb.update_ticks()
cb.set_label('E [V/m]')
pl.title('Singlegrid electric field')
sp1.ticklabel_format(style='sci', scilimits=(0,0),axis='x') 
sp1.ticklabel_format(style='sci', scilimits=(0,0),axis='y')

sp2 = pl.subplot(2,2,2, sharex=sp1, sharey=sp1)
pl.pcolormesh(x_grid_probes, y_grid_probes, 
	np.log10(np.sqrt(Ex_multigrid_matrix**2+Ey_multigrid_matrix**2).T), vmin=vmin, vmax=vmax)
for ii in xrange(pic_multigrid.n_grids):
    sp2.plot(pic_multigrid.pic_list[ii].pic_internal.chamb.Vx, pic_multigrid.pic_list[ii].pic_internal.chamb.Vy, '.-')
pl.xlabel('x [m]')
pl.ylabel('y [m]')
cb=pl.colorbar(); pl.axis('equal')
cb.formatter.set_powerlimits((0, 0))
cb.update_ticks()
cb.set_label('E [V/m]')
pl.title('Multigrid electric field')
sp1.ticklabel_format(style='sci', scilimits=(0,0),axis='x') 
sp1.ticklabel_format(style='sci', scilimits=(0,0),axis='y')

sp3 = pl.subplot(2,2,3, sharex=sp1, sharey=sp1)
pl.pcolormesh(x_grid_probes, y_grid_probes, rho_singlegrid_matrix.T)
for ii in xrange(pic_multigrid.n_grids):
    sp2.plot(pic_multigrid.pic_list[ii].pic_internal.chamb.Vx, pic_multigrid.pic_list[ii].pic_internal.chamb.Vy, '.-')
pl.xlabel('x [m]')
pl.ylabel('y [m]')
cb=pl.colorbar(); pl.axis('equal')
cb.formatter.set_powerlimits((0, 0))
cb.update_ticks()
cb.set_label('Rho')
pl.title('Singlegrid distribution')
sp1.ticklabel_format(style='sci', scilimits=(0,0),axis='x') 
sp1.ticklabel_format(style='sci', scilimits=(0,0),axis='y')

sp4 = pl.subplot(2,2,4, sharex=sp1, sharey=sp1)
pl.pcolormesh(x_grid_probes, y_grid_probes, rho_multigrid_matrix.T)
for ii in xrange(pic_multigrid.n_grids):
    sp2.plot(pic_multigrid.pic_list[ii].pic_internal.chamb.Vx, pic_multigrid.pic_list[ii].pic_internal.chamb.Vy, '.-')
pl.xlabel('x [m]')
pl.ylabel('y [m]')
cb=pl.colorbar(); pl.axis('equal')
cb.formatter.set_powerlimits((0, 0))
cb.update_ticks()
cb.set_label('Rho')
pl.title('Multigrid distribution')
sp1.ticklabel_format(style='sci', scilimits=(0,0),axis='x') 
sp1.ticklabel_format(style='sci', scilimits=(0,0),axis='y')




pl.show()

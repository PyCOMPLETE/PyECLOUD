import sys
sys.path.append('../../../')

import os
import numpy as np
import time

import PyPARIS_sim_class.Simulation as sim_mod

n_pinches_to_average = 3
save_efields = True
save_rho = True
transpose_to_natural_ordering_of_xyz = False

def one_pinch(save_rho=True, save_sigmas=False, save_coords=False):
    
    if os.path.exists('simulation_status.sta'):
        os.remove('simulation_status.sta')
    
    ring = sim_mod.get_serial_CPUring()
    
    list_slice_objects = ring.pieces_to_be_treated # Head is the last element
    list_machine_elements = ring.sim_content.mypart
    # Truck the beam to the first ecloud (excluded)
    for ee in list_machine_elements:
        if ee in ring.sim_content.my_list_eclouds:
            first_ecloud = ee
            break
        for ss in list_slice_objects[::-1]:
            ee.track(ss)
    
    # Record pinch info
    first_ecloud.save_ele_distributions_last_track = True
    first_ecloud.save_ele_field = True
    first_ecloud.save_ele_potential = True
    
    first_ecloud._reinitialize() # Needed to prepare storage space
    
    N_slices = len(list_slice_objects)
    z_centers = []
    t_start = time.mktime(time.localtime())
    for i_ss, ss in enumerate(list_slice_objects[::-1]):
        if np.mod(i_ss, 20)==0:
            print(("%d / %d"%(i_ss, N_slices)))
        first_ecloud.track(ss)
        if ss.slice_info != 'unsliced':
            z_centers.append(ss.slice_info['z_bin_center'])
    
    first_ecloud._finalize()
    z_centers = z_centers[::-1] # HEADTAIL convention
    
    t_end = time.mktime(time.localtime())
    
    print(('Track time %.2f s' % (t_end - t_start)))
    
    dd = {}
    dd['phi'] = first_ecloud.phi_ele_last_track
    if save_rho:
        dd['rho'] = first_ecloud.rho_ele_last_track
    if save_sigmas:
        dd['sigma_x_beam'] = ring.sim_content.bunch.sigma_x()
        dd['sigma_y_beam'] = ring.sim_content.bunch.sigma_y()
        dd['sigma_z_beam'] = ring.sim_content.bunch.sigma_z()
    if save_coords:
        dd['xg'] = first_ecloud.spacech_ele.xg
        dd['yg'] = first_ecloud.spacech_ele.yg
        dd['zg'] = z_centers

    return dd

dd = one_pinch(save_rho=save_rho, save_sigmas=True, save_coords=True)
dd['phi'] /= 1.*n_pinches_to_average
if save_rho:
    dd['rho'] /= 1.*n_pinches_to_average

for i in range(n_pinches_to_average-1):
    print(('Running pinch #%d/%d.'%(i,n_pinches_to_average)))
    dd_temp = one_pinch(save_rho=save_rho, save_sigmas=False, save_coords=False)
    dd['phi'] += dd_temp['phi']/(1.*n_pinches_to_average)
    if save_rho:
        dd['rho'] += dd_temp['rho']/(1.*n_pinches_to_average)

if save_efields:
    dx = dd['xg'][1] - dd['xg'][0]
    dy = dd['yg'][1] - dd['yg'][0]
    dz = dd['zg'][1] - dd['zg'][0]
    dd['Ex'] = np.zeros_like(dd['phi'])
    dd['Ey'] = np.zeros_like(dd['phi'])
    dd['Ez'] = np.zeros_like(dd['phi'])
    
    dd['Ex'][:,1:-1,:] = -0.5/dx*( dd['phi'][:,2:,:] - dd['phi'][:,0:-2,:])
    dd['Ey'][:,:,1:-1] = -0.5/dy*( dd['phi'][:,:,2:] - dd['phi'][:,:,0:-2])
    dd['Ez'][1:-1,:,:] = -0.5/dz*( dd['phi'][2:,::] - dd['phi'][0:-2,:,:])

if transpose_to_natural_ordering_of_xyz:
    dd['phi'] = dd['phi'].transpose(1,2,0)
    dd['rho'] = dd['rho'].transpose(1,2,0)
    dd['Ex']  = dd['Ex'].transpose(1,2,0)
    dd['Ey']  = dd['Ey'].transpose(1,2,0)
    dd['Ey']  = dd['Ey'].transpose(1,2,0)

# Save pinch to file
import scipy.io as sio
sio.savemat('pinch_pic_data.mat', dd, oned_as = 'row')

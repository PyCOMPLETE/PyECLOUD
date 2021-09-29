import os
import numpy as np
import scipy.io as sio
import time

from PyPARIS_sim_class import Simulation as sim_mod
import PyPARIS.util as pu


# start-settings-section
cos_amplitude = 1.00000000e-04
sin_amplitude = 0.00000000e+00
N_oscillations = 3.00000000e+00

flag_no_bunch_charge = False
flag_plots = True

plane = 'y'

sim_param_file = '../reference_simulation/Simulation_parameters.py'
sim_param_amend_files = [
        '../Simulation_parameters_amend.py',
        'Simulation_parameters_amend_for_sin_response.py']
# end-settings-section

# Instantiate simulation
sim_content = sim_mod.Simulation(param_file=sim_param_file)

# Here sim_content.pp can be edited (directly and through files)
for ff in sim_param_amend_files:
    sim_content.pp.update(param_file=ff)

# Add ring of CPU information (mimicking the master core)
pu.get_sim_instance(sim_content,
        N_cores_pretend=sim_content.pp.n_segments,
        id_pretend=sim_content.pp.n_segments-1,
        init_sim_objects_auto=False)
assert(sim_content.ring_of_CPUs.I_am_the_master)

# Initialize machine elements
sim_content.init_all()

# Initialize master to get the beam
if os.path.exists('simulation_status.sta'):
    os.remove('simulation_status.sta')
slices = sim_content.init_master()
N_slices = len(slices)

# Re-center all slices
for ss in slices:
    if ss.macroparticlenumber:
        ss.x -= ss.mean_x()
        ss.xp -= ss.mean_xp()
        ss.y -= ss.mean_y()
        ss.yp -= ss.mean_yp()

# Optionally remove charge from bunch
if flag_no_bunch_charge:
    for ss in slices:
        ss.particlenumber_per_mp = 1e-10

# Get slice centers
z_slices = np.array([ss.slice_info['z_bin_center'] for ss in slices])

# Get z_step beween slices and define z_range
z_step = z_slices[1] - z_slices[0]
z_range = z_slices[-1] - z_slices[0] + z_step # Last term is to make the sampled 
                                              # sinusoids more orthogonal
# Generate ideal sinusoidal distortion 
r_ideal = (sin_amplitude * np.sin(2*np.pi*N_oscillations*z_slices/z_range)
         + cos_amplitude * np.cos(2*np.pi*N_oscillations*z_slices/z_range))

# Add sinusoidal distortion to particles
for ss in slices:
    if ss.macroparticlenumber>0:
        #if ss.mean_z() < 0:
        if plane == 'x':
            ss.x += sin_amplitude * np.sin(2*np.pi*N_oscillations*ss.z/z_range)
            ss.x += cos_amplitude * np.cos(2*np.pi*N_oscillations*ss.z/z_range)
        elif plane == 'y':
            ss.y += sin_amplitude * np.sin(2*np.pi*N_oscillations*ss.z/z_range)
            ss.y += cos_amplitude * np.cos(2*np.pi*N_oscillations*ss.z/z_range)
        else:
            raise ValueError('What?!')

# Measure
if plane == 'x':
    r_slices = np.array([ss.mean_x() for ss in slices])
elif plane == 'y':
    r_slices = np.array([ss.mean_y() for ss in slices])
int_slices = np.array([ss.intensity for ss in slices])

# Simulate e-cloud interactions
t_start = time.mktime(time.localtime())
dpr_slices = []
rho_slices = []
for i_ss, ss in enumerate(slices[::-1]):
    if np.mod(i_ss, 20)==0:
        print(("%d / %d"%(i_ss, N_slices)))
    for i_ee, ee in enumerate(sim_content.parent_eclouds):
        ee.track(ss)
        if i_ee == 0:
            temp_rho = ee.cloudsim.cloud_list[0].rho.copy()
        else:
            temp_rho += ee.cloudsim.cloud_list[0].rho.copy()
    if plane == 'x':
        dpr_slices.append(ss.mean_xp())
    elif plane == 'y':
        dpr_slices.append(ss.mean_yp())
    rho_slices.append(temp_rho)
dpr_slices = np.array(dpr_slices[::-1])
rho_slices = np.array(rho_slices[::-1])
t_end = time.mktime(time.localtime())
print(('Ecloud sim time %.2f s' % (t_end - t_start)))

dpr_slices_all_clouds = dpr_slices * sim_content.n_segments

# Savings and plots
first_ecloud = sim_content.parent_eclouds[0]
xg = first_ecloud.cloudsim.spacech_ele.xg
yg = first_ecloud.cloudsim.spacech_ele.yg

if plane == 'x':
    i_yzero = np.argmin(np.abs(yg))
    rho_cut = rho_slices[:, :, i_yzero]
    rg = xg
elif plane == 'y':
    i_xzero = np.argmin(np.abs(xg))
    rho_cut = rho_slices[:, i_xzero, :]
    rg = yg


sio.savemat('response.mat',{
    'plane': plane,
    'z_slices': z_slices,
    'r_slices': r_slices,
    'r_ideal': r_ideal,
    'dpr_slices_all_clouds': dpr_slices_all_clouds,
    'xg': xg,
    'yg': yg,
    'rho_cut': rho_cut,
    })

if flag_plots:
    import matplotlib.pyplot as plt
    plt.close('all')
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(3,1,1)
    ax2 = fig1.add_subplot(3,1,2, sharex=ax1)
    ax3 = fig1.add_subplot(3,1,3, sharex=ax1)

    ax1.plot(z_slices, int_slices)
    ax2.plot(z_slices, r_slices)
    ax3.plot(z_slices, dpr_slices)

    for ax in [ax1, ax2, ax3]:
        ax.grid(True)


    fig2 = plt.figure(2)
    ax21 = fig2.add_subplot(2,1,1)
    ax22 = fig2.add_subplot(2,1,2, sharex=ax21)

    ax21.pcolormesh(z_slices, rg, rho_cut.T)
    ax21.plot(z_slices, r_slices, 'k', lw=2)
    ax22.plot(z_slices, dpr_slices)
    ax22.set_ylim(np.nanmax(np.abs(dpr_slices))*np.array([-1, 1]))
    ax22.grid(True)
    plt.show()




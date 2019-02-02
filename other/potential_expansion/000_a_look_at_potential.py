import sys, os
sys.path.append(os.path.expanduser('../../../'))
sys.path.append(os.path.expanduser('../../../PyHEADTAIL/'))

from scipy.constants import c, e, m_p
import numpy as np
import pylab as plt
import PyECLOUD.myfilemanager as mfm
import PyECLOUD.mystyle as ms

n_segments = 1
machine_configuration = '6.5_TeV_collision_tunes'
p0_GeV = 2000

z_cut = 2.5e-9 * c
n_slices = 150
L_ecloud = 1000.

sigma_z = 10e-2
epsn_x = 2.5e-6
epsn_y = 2.5e-6

sparse_solver = 'PyKLU'  # 'scipy_slu'


# Define the machine
#============================
from machines_for_testing import LHC
machine = LHC(machine_configuration=machine_configuration,
              optics_mode='smooth', n_segments=n_segments, p0=p0_GeV * 1e9 * e / c)

bunch = machine.generate_6D_Gaussian_bunch(
    n_macroparticles=3000000, intensity=1e11,
    epsn_x=epsn_x, epsn_y=epsn_y, sigma_z=sigma_z)


from PyHEADTAIL.particles.slicing import UniformBinSlicer
slicer = UniformBinSlicer(n_slices=n_slices, z_cuts=(-z_cut, z_cut) )
import PyECLOUD.PyEC4PyHT as PyEC4PyHT

ecloud_multigrid = PyEC4PyHT.Ecloud(
    L_ecloud=L_ecloud, slicer=slicer,
    Dt_ref=20e-12, pyecl_input_folder='./pyecloud_config_LHC',
    chamb_type='polyg' ,
    filename_chm='LHC_chm_ver.mat', Dh_sc=1e-3,
    init_unif_edens_flag=1,
    init_unif_edens=1e7,
    N_mp_max=3000000,
    nel_mp_ref_0 = 1e7*np.pi*2e-2**2 / (0.7 * 3000000),
    B_multip=[0.],
    PyPICmode='ShortleyWeller_WithTelescopicGrids',
    f_telescope=0.3,
    target_grid={'x_min_target': -5 * bunch.sigma_x(), 'x_max_target': 5 * bunch.sigma_x(),
                 'y_min_target': -5 * bunch.sigma_y(), 'y_max_target': 5 * bunch.sigma_y(),
                 'Dh_target': .2 * bunch.sigma_x()},
    N_nodes_discard=10.,
    N_min_Dh_main=10,
    sparse_solver=sparse_solver)



ecloud_multigrid.save_ele_distributions_last_track = True
ecloud_multigrid.save_ele_field = True
ecloud_multigrid.save_ele_potential = True

ecloud_multigrid.track(bunch)

slices = bunch.get_slices(ecloud_multigrid.slicer)

xg = ecloud_multigrid.spacech_ele.xg
yg = ecloud_multigrid.spacech_ele.yg

import scipy.io as sio

sio.savemat('pinch_pic_data.mat', {
    'xg': xg,
    'yg': yg,
    'zg': slices.z_centers,
    'rho': ecloud_multigrid.rho_ele_last_track,
    'phi': ecloud_multigrid.phi_ele_last_track,
    'Ex': ecloud_multigrid.Ex_ele_last_track,
    'Ey': ecloud_multigrid.Ey_ele_last_track,
    }, oned_as = 'row')
    

x_obs = 0. 
ix_obs = np.argmin(np.abs(xg - x_obs))


plt.close('all')
ms.mystyle_arial()
fig1 = plt.figure()
ax1 = fig1.add_subplot(3,1,1)
ax2 = fig1.add_subplot(3,1,2, sharex=ax1)
ax3 = fig1.add_subplot(3,1,3, sharex=ax1)

ax1.pcolormesh(slices.z_centers, yg, 
        ecloud_multigrid.rho_ele_last_track[:,ix_obs,:].T)
ax2.pcolormesh(slices.z_centers, yg, 
        ecloud_multigrid.Ey_ele_last_track[:,ix_obs,:].T)
ax3.pcolormesh(slices.z_centers, yg, 
        ecloud_multigrid.phi_ele_last_track[:,ix_obs,:].T)



plt.show()


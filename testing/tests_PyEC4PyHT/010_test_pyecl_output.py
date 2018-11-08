import sys, os
sys.path.append(os.path.expanduser('../../../'))
sys.path.append(os.path.expanduser('../../../PyHEADTAIL/'))

from scipy.constants import c, e, m_p
import numpy as np
import pylab as pl
import myfilemanager as mlm
import PyECLOUD.mystyle as ms

n_segments =1
machine_configuration = '6.5_TeV_collision_tunes'
p0_GeV = 2000

z_cut = 2.5e-9*c
n_slices = 300
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


import PyECLOUD.PyEC4PyHT as PyEC4PyHT


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
        sparse_solver = sparse_solver,
        save_pyecl_outp_as = 'test_saving',
        Dt = 25e-12,#needed for saving
        )

print 'Track_bunch'
ecloud_multigrid.track(bunch)
print 'Done.'

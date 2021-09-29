import sys, os
sys.path.append(os.path.expanduser('../../../'))
sys.path.append(os.path.expanduser('../../../PyHEADTAIL/'))

import numpy as np

from PyHEADTAIL.particles.slicing import UniformBinSlicer

n_segments = 10
N_turns = 1

epsn_x = 2.5e-6
epsn_y = 2.5e-6

init_unif_edens_flag = 1
init_unif_edens = 1e11
N_MP_ele_init = 100000
N_mp_max = N_MP_ele_init * 4.

# define the machine
from machines_for_testing import SPS
machine = SPS(n_segments=n_segments, machine_configuration='Q26-injection')

# compute sigma x and y
inj_optics = machine.transverse_map.get_injection_optics()
sigma_x = np.sqrt(inj_optics['beta_x'] * epsn_x / machine.betagamma)
sigma_y = np.sqrt(inj_optics['beta_y'] * epsn_y / machine.betagamma)

# define apertures and Dh_sc to simulate headtail conditions
x_aper = 20 * sigma_x
y_aper = 20 * sigma_y
Dh_sc = 2 * x_aper / 128

# define MP size
nel_mp_ref_0 = init_unif_edens * 4 * x_aper * y_aper / N_MP_ele_init

# define an electron cloud
import PyECLOUD.PyEC4PyHT as PyEC4PyHT
from PyHEADTAIL.particles.slicing import UniformBinSlicer
slicer = UniformBinSlicer(n_slices=64, n_sigma_z=2.)
ecloud = PyEC4PyHT.Ecloud(L_ecloud=machine.circumference / n_segments, slicer=slicer ,
                          Dt_ref=25e-12, pyecl_input_folder='./drift_sim',
                          x_aper=x_aper, y_aper=y_aper, Dh_sc=Dh_sc,
                          init_unif_edens_flag=init_unif_edens_flag,
                          init_unif_edens=init_unif_edens,
                          N_mp_max=N_mp_max,
                          nel_mp_ref_0=nel_mp_ref_0)


# install ecloud in the machine
machine.install_after_each_transverse_segment(ecloud)

# setup transverse losses (to "protect" the ecloud)
import PyHEADTAIL.aperture.aperture as aperture
apt_xy = aperture.EllipticalApertureXY(x_aper=ecloud.cloudsim.chamb.x_aper, y_aper=ecloud.cloudsim.chamb.y_aper)
machine.one_turn_map.append(apt_xy)

# generate a bunch
bunch = machine.generate_6D_Gaussian_bunch(n_macroparticles=300000, intensity=1.5e11, epsn_x=epsn_x, epsn_y=epsn_y, sigma_z=.11)

# simulate
for i_turn in range(N_turns):
	print('Turn', i_turn)
	machine.track(bunch, verbose=True)



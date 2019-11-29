import numpy as np
from scipy.constants import e as qe

from PyPARIS_sim_class import LHC_custom
from PyECLOUD import PyEC4PyHT
from PyHEADTAIL.particles.slicing import UniformBinSlicer

import Simulation_parameters as pp

n_mp_bunch = 100000
intensity = 1.2e11
epsn_x = 2.5e-6
epsn_y = 2.5e-6
sigma_z = 9e-2

#################
# Build machine #
#################

machine = LHC_custom.LHC(n_segments=1,
        machine_configuration='HLLHC-injection'
        )

################
# Make a bunch #
################

bunch = machine.generate_6D_Gaussian_bunch(
        n_mp_bunch, intensity, epsn_x, epsn_y, sigma_z)


##########
# Slicer #
##########

slicer = UniformBinSlicer(
        n_slices = pp.n_slices, z_cuts=(-pp.z_cut, pp.z_cut))

#########
# Slice #
#########

slices = bunch.extract_slices(slicer)


##################
# Make an ecloud #
##################

inj_opt = machine.transverse_map.get_injection_optics()
beta_x_smooth = inj_opt['beta_x']
beta_y_smooth = inj_opt['beta_y']
sigma_x_smooth = np.sqrt(beta_x_smooth*epsn_x/machine.betagamma)
sigma_y_smooth = np.sqrt(beta_y_smooth*epsn_y/machine.betagamma)

target_grid_arcs = {
    'x_min_target':-pp.target_size_internal_grid_sigma*sigma_x_smooth,
    'x_max_target':pp.target_size_internal_grid_sigma*sigma_x_smooth,
    'y_min_target':-pp.target_size_internal_grid_sigma*sigma_y_smooth,
    'y_max_target':pp.target_size_internal_grid_sigma*sigma_y_smooth,
    'Dh_target':pp.target_Dh_internal_grid_sigma*sigma_x_smooth}

nel_mp_ref_0 = pp.init_unif_edens_dip*4*pp.x_aper*pp.y_aper/pp.N_MP_ele_init_dip

ecloud = PyEC4PyHT.Ecloud(slice_by_slice_mode=True,
    L_ecloud=1., slicer=None,
    Dt_ref=pp.Dt_ref, pyecl_input_folder=pp.pyecl_input_folder,
    chamb_type = pp.chamb_type,
    x_aper=pp.x_aper, y_aper=pp.y_aper,
    filename_chm=pp.filename_chm,
    PyPICmode = pp.PyPICmode,
    Dh_sc=pp.Dh_sc_ext,
    N_min_Dh_main = pp.N_min_Dh_main,
    f_telescope = pp.f_telescope,
    N_nodes_discard = pp.N_nodes_discard,
    target_grid = target_grid_arcs,
    init_unif_edens_flag=pp.init_unif_edens_flag_dip,
    init_unif_edens=pp.init_unif_edens_dip,
    N_mp_max=pp.N_mp_max_dip,
    nel_mp_ref_0=nel_mp_ref_0,
    B_multip=pp.B_multip_dip,
    enable_kick_x = pp.enable_kick_x,
    enable_kick_y = pp.enable_kick_y)



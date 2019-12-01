import numpy as np
from scipy.constants import e as qe
import matplotlib.pyplot as plt
from PyPARIS_sim_class import LHC_custom
from PyECLOUD import PyEC4PyHT
from PyHEADTAIL.particles.slicing import UniformBinSlicer

import Simulation_parameters as pp

plt.close('all')

n_mp_slice = 25000
intensity = 1.2e11
epsn_x = 2.5e-6
epsn_y = 2.5e-6
sigma_z = 9.7e-2
n_slices = 150
n_mp_bunch = n_mp_slice*n_slices
long_unif = True

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
if long_unif:
    bunch.z = 4*sigma_z*(np.random.rand(len(bunch.z))-0.5)

##########
# Slicer #
##########

slicer = UniformBinSlicer(
        n_slices = n_slices, z_cuts=(-pp.z_cut, pp.z_cut))

#########
# Slice #
#########

slices = bunch.extract_slices(slicer)
z_slices = np.array([
    ss.slice_info['z_bin_center'] for ss in slices])
intensity_slices = np.array([
    ss.intensity for ss in slices])


fig1 = plt.figure(1)
axl = fig1.add_subplot(111)
axl.plot(z_slices, intensity_slices)

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
#    flag_reinterp_fields_at_substeps=False,
    force_interp_at_substeps_interacting_slices=False,
    #Dt_ref=1.,
    Dt_ref=pp.Dt_ref,
    pyecl_input_folder=pp.pyecl_input_folder,
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
    B_multip=[0.],
    enable_kick_x = pp.enable_kick_x,
    enable_kick_y = pp.enable_kick_y,
    kick_mode_for_beam_field=False)

sc = ecloud.beam_PyPIC_state.pic_internal
i_y0 = np.argmin(np.abs(sc.yg))


mpe = ecloud.cloudsim.cloud_list[0].MP_e
mpe.x_mp[0] = .1e-3 
mpe.vx_mp[0] = 0.
mpe.y_mp[0] = 0.
mpe.vy_mp[0] = 0.
mpe.z_mp[0] = 0.
N_track = 1
x_record = []
field_record = []
for ss in slices[::-1]:
    #ecloud.track(slices[int(n_slices/2)])
    ecloud.track(ss)
    field_record.append(sc.efx[:, i_y0].copy())
    x_record.append(mpe.x_mp[:N_track].copy())
x_record = np.array(x_record[::-1])

fig2 = plt.figure(2)
axx = fig2.add_subplot(111)
axx.plot(z_slices, x_record, '.-')
axx.grid(True)
plt.show()


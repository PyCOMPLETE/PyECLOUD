"""Track one bunch with multigrid and record the electron cloud on the finest grid.

Based on 008_test_multigrid_pinch.py. Requires installed PyECLOUD, PyPIC and
PyHEADTAIL. Input paths are relative to this script. Diagnostics stay in memory;
interactive plots show electron number density in the x=0, y=0 and z=0 planes.

Ecloud's built-in diagnostics retain the innermost grid, not every refinement
level. Arrays have shape (n_slices, len(x_grid), len(y_grid)), in increasing z
order (the bunch is tracked from positive to negative z). A subsequent track()
replaces Ecloud's logs. The recorded arrays contain the electron contribution only.
"""

from pathlib import Path
from time import perf_counter

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c, e
from scipy.interpolate import interp1d
from scipy.io import loadmat

from PyHEADTAIL.particles.slicing import UniformBinSlicer
from PyECLOUD.PyEC4PyHT import Ecloud
from machines_for_testing import LHC


input_dir = Path(__file__).resolve().parent
np.random.seed(12345)

# Same physical parameters as 008, with more slices and electron macroparticles.
n_macroparticles = 3_000_000
n_slices = 301
p0_GeV = 2000
z_cut = 2.5e-9 * c
L_ecloud = 1000.
sigma_z = 0.10
epsn_x = epsn_y = 2.5e-6
init_unif_edens = 1e7
N_electron_macroparticles = 1_000_000  # target initial count inside the chamber
N_mp_max = 3_000_000

# The 2D cloud uses macroparticle weights in electrons per metre. Include the
# polygon area when converting the physical volume density to a target count.
# Initialization samples a bounding rectangle and rejects points outside the
# chamber, so the actual accepted count fluctuates slightly around the target.
chamber_file = input_dir / 'LHC_chm_ver.mat'
chamber_data = loadmat(chamber_file)
vx = chamber_data['Vx'].ravel()
vy = chamber_data['Vy'].ravel()
chamber_area = 0.5 * abs(np.dot(vx, np.roll(vy, -1))
                         - np.dot(vy, np.roll(vx, -1)))
nel_mp_ref_0 = init_unif_edens * chamber_area / N_electron_macroparticles

machine = LHC(
    machine_configuration='6.5_TeV_collision_tunes',
    optics_mode='smooth', n_segments=1, p0=p0_GeV * 1e9 * e / c,
)
bunch = machine.generate_6D_Gaussian_bunch(
    n_macroparticles=n_macroparticles, intensity=1e11,
    epsn_x=epsn_x, epsn_y=epsn_y, sigma_z=sigma_z,
)
sigma_x = bunch.sigma_x()
sigma_y = bunch.sigma_y()
slicer = UniformBinSlicer(n_slices=n_slices, z_cuts=(-z_cut, z_cut))

ecloud = Ecloud(
    L_ecloud=L_ecloud, slicer=slicer, Dt_ref=20e-12,
    pyecl_input_folder=str(input_dir / 'pyecloud_config_LHC'),
    chamb_type='polyg', filename_chm=str(chamber_file),
    init_unif_edens_flag=1, init_unif_edens=init_unif_edens,
    N_mp_max=N_mp_max,
    nel_mp_ref_0=nel_mp_ref_0,
    B_multip=[0.], x_beam_offset=0., y_beam_offset=0.,
    sparse_solver='scipy_slu',
    PyPICmode='ShortleyWeller_WithTelescopicGrids',
    Dh_sc=1e-3, f_telescope=0.3,
    target_grid={
        'x_min_target': -5 * sigma_x, 'x_max_target': 5 * sigma_x,
        'y_min_target': -5 * sigma_y, 'y_max_target': 5 * sigma_y,
        'Dh_target': 0.2 * sigma_x,
    },
    N_nodes_discard=10, N_min_Dh_main=10,
)
print(f'Initial electron macroparticles: {ecloud.cloudsim.cloud_list[0].MP_e.N_mp:,}'
      f' (target {N_electron_macroparticles:,})')

# Set these after construction. track() calls _reinitialize() itself to prepare
# the storage, then _finalize() converts the snapshots to arrays ordered by z.
ecloud.save_ele_distributions_last_track = True
ecloud.save_ele_potential_and_field = True

t_start = perf_counter()
ecloud.track(bunch)
print(f'Multigrid tracking time: {perf_counter() - t_start:.3f} s')

# Keep convenient aliases for interactive use, with explicit physical units.
x_grid = ecloud.spacech_ele.xg.copy()  # m; finest grid, including its boundary
y_grid = ecloud.spacech_ele.yg.copy()  # m
z_centers = bunch.get_slices(slicer).z_centers.copy()  # m, increasing order
rho_ele = ecloud.rho_ele_last_track  # C/m^3 (negative for electrons)
n_ele = -rho_ele / e  # electron number density, m^-3
phi_ele = ecloud.phi_ele_last_track  # V
Ex_ele = ecloud.Ex_ele_last_track  # V/m
Ey_ele = ecloud.Ey_ele_last_track  # V/m

print(f'Recorded grid arrays with shape {rho_ele.shape} (slice, x, y)')

# Interpolate onto the exact zero planes (also works with an even slice count).
# The z coordinate labels successive cloud snapshots during
# the bunch passage, not a simultaneous 3D cloud distribution.
n_at_x0 = interp1d(x_grid, n_ele, axis=1)(0.)  # (z, y)
n_at_y0 = interp1d(y_grid, n_ele, axis=2)(0.)  # (z, x)
n_at_z0 = interp1d(z_centers, n_ele, axis=0)(0.)  # (x, y)

plt.ion()
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), layout='constrained')
vmax = max(n_at_x0.max(), n_at_y0.max(), n_at_z0.max())
for ax, horizontal, vertical, density, xlabel, ylabel, title in (
    (axes[0], z_centers * 1e2, y_grid * 1e3, n_at_x0.T,
     'z [cm]', 'y [mm]', 'x = 0'),
    (axes[1], z_centers * 1e2, x_grid * 1e3, n_at_y0.T,
     'z [cm]', 'x [mm]', 'y = 0'),
    (axes[2], x_grid * 1e3, y_grid * 1e3, n_at_z0.T,
     'x [mm]', 'y [mm]', 'z = 0'),
):
    mesh = ax.pcolormesh(horizontal, vertical, density, shading='auto',
                         vmin=0., vmax=vmax, cmap='viridis')
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
axes[2].set_aspect('equal')
fig.colorbar(mesh, ax=axes, label=r'Electron number density [m$^{-3}$]')
fig.suptitle('Electron cloud pinch — finest grid')
# Keep the GUI open when launched as a script; its toolbar supports zoom/pan.
plt.show(block=True)

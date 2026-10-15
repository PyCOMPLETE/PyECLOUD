# Single-bunch simulation parameters

The file Simulation_parameters.py allows configuring your PyECLOUD-PyHEADTAIL single-bunch simulation.
In the following we describe its different parts:

## Initial imports
In the beginning of the file you can import python modules or objects required to configure the simulation. In this example we import some constants:

```python
from scipy.constants import c
from scipy.constants import m_p
from scipy.constants import e as qe
```

## Description of the accelerator
By choosing ```"Synchrotron"``` as machine class it is possible to define all parameters of the ring (this is the recommended way of configuring the simulation):
```python
machine_class = 'Synchrotron'
```
The description of the accelerator is defined by providing the following parameters that are required to define a PyHEADTAIL synchrotron. See the [PyHEADTAIL Synchrotron implementation and constructor documentation](https://github.com/PyCOMPLETE/PyHEADTAIL/blob/master/PyHEADTAIL/machines/synchrotron.py) for their definitions.
```python
optics_mode = 'smooth'
charge = qe
mass = m_p
p0 = 450e9 * qe / c
circumference = 26658.8832
n_segments = 16
name = None
s = None
alpha_x = 0.
beta_x = 92.7
D_x = 0.
alpha_y = 0
beta_y = 93.2
D_y = 0.
accQ_x = 62.27
accQ_y = 60.295
Qp_x = 0.
Qp_y = 0.
app_x = 0.
app_y = 0.
app_xy = 0.
longitudinal_mode = 'non-linear'
Q_s = None
alpha_mom_compaction = 3.225e-04
h_RF = 35640
V_RF = 5e6
dphi_RF = 0.
p_increment = 0.
RF_at = 'end_of_transverse' # needs to be at the end to be compatible with PyPARIS parallelization
wrap_z = False
other_detuners = []
```

The parameter ```n_non_parallelizable``` (used in previous versions) is now handled automatically and should not be added to the input file.

## Transverse feedback
An ideal bunch-by-bunch transverse feedback can be enabled. The damping rates are defined as the time constant observed on the envelope of the bunch centroid motion in the presence of the feedback alone:
```python
# Transverse Damper Settings
enable_transverse_damper = False
dampingrate_x = 100.
dampingrate_y = 100.
```

## Beam parameters

The parameters in this section define the bunch.

Optionally the bunch (particle-by-particle coordinates) can be loaded from file using the following input variable:
```python
bunch_from_file = None
```
An example of how to save and load a bunch can be found in `PyPARIS_sim_class/examples`.

Alternatively, the bunch can be defined using the following input
```python
intensity = 1.2e+11    # Number of particles in the bunch
epsn_x = 2.5e-6        # Normalized horizontal r.m.s emittance [m]
epsn_y = 2.5e-6        # Normalized vertical r.m.s emittance [m]
sigma_z = 9.181144e-02 # R.m.s bunch length [m]
```
The generated bunch is matched in all planes.

An initial transverse displacement can be defined by the following parameters:
```python
x_kick_in_sigmas = 0.1
y_kick_in_sigmas = 0.1
```
The average position and angle of the bunch longitudinal slices can be shifted to zero (to correct for random fluctuations) by setting:
```python
recenter_all_slices = True
```
## Longitudinal slicing
The beam is sliced longitudinally to compute its interaction with the e-cloud:
```python
n_slices = 200
z_cut = 2.5e-9/2*c # For slicing
```
```z_cut``` defines the portion of the bunch that is sliced, which is [-z_cut, z_cut].

## Number of macroparticles in the bunch
The number of macroparticles in the bunch is defined by the following parameters:
```
macroparticles_per_slice = 5000
n_macroparticles = macroparticles_per_slice*n_slices
```
## Multijob setup
For very long simulations it is convenient to split the simulations over several shorter jobs. This also allows recovering the simulation from the end of the last successful job in case of problems. This is controlled through the following parameters:

```python
N_turns = 128 # Per job
N_turns_target = 20000 # Entire simulation
```  
For clusters supporting job resubmissions the job can be resubmitted automatically by setting:
```python
check_for_resubmit = True
```
In single-bunch multijob mode, `check_for_resubmit` defaults to `True` when omitted.
Set it explicitly to `False` to manage job resubmission yourself.
Currently supported submission systems are:

- HTCondor
- SLURM

Select the submission system using the optional parameter (names are case-insensitive):
```python
submission_system = 'HTCondor' # Default; alternatively 'SLURM'
```
For `HTCondor`, the simulation exits with code 177. See the
[Batch Jobs guide on DAGMan continuation](../batch/stability.md#dagman-continuation)
for how to handle this code and submit the next simulation part.
For `Slurm`, it runs `sbatch job.cmd` from the simulation's working directory;
that directory must contain your submission script under that name.
To run a custom command instead, configure:
```python
resubmit_command = 'sbatch my_job.sh'
```
This parameter is set to `None` by default and if configured, it takes precedence over `submission_system`.
The custom command is interpreted by a shell. The simulation process exits with
the command's return code. PyPARIS's multiprocessing launcher preserves code 177
but reports other worker failures as code 1.

The information about the multijob run is kept in the file ```simulation_status.sta```. This file needs to be removed to start the simulation from scratch.

## Stop criteria
The simulation can be ended in case a certain fraction of the initial intensity is lost:

```python
sim_stop_frac = 0.9 # Simulation stopped if 10% of the initial intensity is lost
```

Transverse emittance blow-up can also be used as a criterion for ending the simulation:

```python
flag_check_emittance_growth = True
epsn_x_max_growth_fraction = 0.5 # Stop on 50% horizontal emittance blow-up
epsn_y_max_growth_fraction = 0.5 # Stop on 50% vertical emittance blow-up
```

## Saving settings

The output data is buffered and dumped to file at regular intervals. The length of such intervals can be specified by setting:

```python
write_buffer_every = 3
slice_stats_to_store = ['mean_x', 'mean_z']
```

## Footprint mode
The simulation mode can be changed to compute the tune footprint in the presence of e-cloud using the ```footprint_mode``` flag. In this case the instability simulation is not performed and only the footprint data is produced. A bunch with a very large number of macroparticles is used to compute the e-cloud field maps. The footprint particle tunes are measured using a smaller number of macroparticles. For the computation of the footprint the longitudinal motion is automatically switched off.

The footprint mode is enabled and controlled using the following parameters:

```python
footprint_mode = False
n_macroparticles_for_footprint_map = 500000
n_macroparticles_for_footprint_track = 5000
```


## Electron cloud settings
Two kinds of e-cloud interaction can be installed at the end of each machine segment to model the e-cloud in dipoles and quadrupoles.

Most of the e-cloud parameters are defined by PyECLOUD input files in a specified folder:
```python
pyecl_input_folder = './pyecloud_config'
```
Whenever a parameter in the input files is also mentioned in the present file, the value specified here takes priority.

The chamber geometry can be redefined by the following parameters (see the [machine parameters](machine.md#chamber-profile) for more details):

```python
chamb_type = 'polyg'
x_aper = 2.300000e-02
y_aper = 1.800000e-02
filename_chm = 'LHC_chm_ver.mat'
```

The target duration of the electron tracking time sub-steps is defined by:
```
Dt_ref = 5.000000e-12
```

The Particle-In-Cell solver used to compute beam and electron fields is configured by the following parameters (more details can be found in the [simulation parameters](simulation.md#multigrid) and in the [physical models and numerical algorithms](../physics.md)):

```python
PyPICmode = 'ShortleyWeller_WithTelescopicGrids'
N_min_Dh_main = 10.
Dh_sc_ext = .8e-3
f_telescope = 0.3
N_nodes_discard = 5.
target_size_internal_grid_sigma = 10.
target_Dh_internal_grid_sigma = 0.2
custom_target_grid_arcs = None

# # Uncomment for custom grid
# custom_target_grid_arcs = {
#     'x_min_target': -3e-3,
#     'x_max_target': 3e-3,
#     'y_min_target': -3.1e-3,
#     'y_max_target': 3.1e-3,
#     'Dh_target': 7e-5}
```

By setting:
```python
force_interp_at_substeps_interacting_slices = True
```
the electric field acting on the electrons is re-interpolated at each tracking sub-step. If this parameter is set to false the field is evaluated only once for each beam-slice.

The forces of the e-cloud on the beam particle in the transverse planes can be enabled/disabled using the following flags:

```python
enable_kick_x = True
enable_kick_y = False
```

### Dedicated settings for e-cloud in the dipole magnets

The following parameters configure the e-cloud in the dipoles (which are simulated starting with electrons at rest uniformly distributed in the chamber):

```python
enable_arc_dip = True               # Activate interaction with e-cloud in the dipoles
fraction_device_dip = 0.65          # Fraction of the machine circumference with e-cloud in the dipoles
init_unif_edens_flag_dip = 1        # Activate initial uniform distribution for electrons
init_unif_edens_dip = 1.0e+12       # Initial electron density (e-/m^3)
N_MP_ele_init_dip = 500000          # Number of macroparticles
N_mp_max_dip = N_MP_ele_init_dip*4  # Size of arrays used to store macroparticle coordinates
B_multip_dip = [0.5] #T             # Magnetic field (in Tesla)
```

### Dedicated settings for e-cloud in the quadrupole magnets
The following parameters configure the e-cloud in the quadrupole magnets (which are simulated starting with an electron distribution loaded from file):

```python
enable_arc_quad = False             # Activate interaction with e-cloud in the quadrupoles
fraction_device_quad = 7.0e-02      # Fraction of the machine circumference with e-cloud in the quadrupoles

N_mp_max_quad = 2000000             # Size of arrays used to store macroparticle coordinates
B_multip_quad = [0., 12.1] #T       # The second element of the list is the magnetic field gradient (in Tesla/m)
```
The following parameters are used to define the file containing the initial e-cloud distribution:
```python
folder_path = '../../LHC_ecloud_distrib_quads/'
N_mp_ele_quad = 500000
sey_load_quad = 1.3
filename_state = 'combined_distribution_sey_%.2f_sigmat_%.3fns_450Gev_N_mp_%d_symm'%(sey_load_quad, sigma_z/c*1e9,N_mp_ele_quad)
filename_init_MP_state_quad = folder_path + filename_state
```

### Expert e-cloud kicks

Individual e-cloud interactions can be provided at defined kick locations when the optics mode is ```'non-smooth'``` (incompatible with the settings above), using the following parameters:

```python
# Dedicated Kick Element Settings
enable_eclouds_at_kick_elements = False
path_buildup_simulations_kick_elements = 'path_to_myfolder/simulations_PyECLOUD/!!!NAME!!!_sey1.35'
name_MP_state_file_kick_elements = 'MP_state_9.mat'
orbit_factor = 6.250000e-01

```

### Impedance

A simple resonator impedance can be included in the simulation by using the following parameters:

```python
enable_impedance = True
resonator_R_shunt = 25e6
resonator_frequency = 2e9
resonator_Q = 1.
n_slices_wake = 500
```

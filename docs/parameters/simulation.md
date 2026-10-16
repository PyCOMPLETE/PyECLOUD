# Simulation parameters

This input file must be named `simulation_parameters.input`.

In the tables below, MP stands for *macroparticle*, and `Nref` is the reference MP size: the
number of electrons each MP represents. Several checks run at every bunch passage, i.e. at
times $t = n \cdot$ `b_spac` with $n$ an integer.

## Other input files

These parameters give the names of the other three input files, which define the physical model
of the simulation.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `machine_param_file`
  - Name of the [machine parameter](machine.md) file.
* - `secondary_emission_parameters_file`
  - Name of the [secondary emission](secondary-emission.md) configuration file.
* - `beam_parameters_file`
  - Name of the [beam parameter](beam.md) file. In a simulation with several beams, this is the
    master beam.
```

## Secondary beams

The code can simulate electron cloud buildup with more than one circulating beam. To do so,
provide a list of secondary beam files, one per beam. The master beam determines the length of
the simulation, the energy used to calculate the bending field, and the bunch spacing used for
regeneration and saving.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `secondary_beams_file_list`
  - *Optional, default `[]`.* List of secondary beam parameter file names.
```

## Log and progress files

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `logfile_path`
  - Path of a text file that summarizes the running simulation: passage number, number of MPs,
    and number of electrons in the chamber.
* - `progress_path`
  - Path of a text file that reports the simulation progress in percent.
```

## Time sampling

The simulated time interval is set by the length of the beam profile (the number of bunch
passages) given in the beam description.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `Dt`
  - \[s\] Simulation time step. Ignored if the beam profile is imported from a file.
* - `t_end`
  - \[s\] Extra time simulated after the end of the beam profile, during which only the coasting
    beam is present. A negative value stops the simulation before the end of the profile.
```

## Negligible beam linear density

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `lam_th`
  - \[beam particles/m\] Below this beam line density, no primary electrons are generated and
    beam forces on the electrons are neglected.
```

## MP management

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `N_mp_max`
  - Size of the arrays that store MP sizes, positions, and velocities. Choose the regeneration
    settings so that this number of MPs is never exceeded.
* - `N_mp_regen`
  - Upper regeneration threshold. When the number of MPs exceeds `N_mp_regen`, the MP set is
    regenerated with a new `Nref` chosen to reach `N_mp_after_regen` MPs. Checked at every bunch
    passage.
* - `N_mp_regen_low`
  - Lower regeneration threshold. When the number of MPs falls below `N_mp_regen_low`, the MP set
    is regenerated with a new `Nref` chosen to reach `N_mp_after_regen` MPs. Checked at every
    bunch passage.
* - `t_ON_regen_low`
  - \[s\] The lower regeneration check is skipped for $t <$ `t_ON_regen_low`.
* - `N_mp_after_regen`
  - Target number of MPs after a regeneration.
* - `fact_split`
  - MP split factor. When secondary emission occurs and the emitted charge is smaller than
    `fact_split * Nref`, the impacting MP is simply rescaled by the SEY. Otherwise new MPs are
    generated to keep MP sizes as close as possible to `Nref`.
* - `fact_clean`
  - MP clean factor. MPs smaller than `fact_clean * Nref` are removed from the MP set. Checked at
    every bunch passage.
* - `nel_mp_ref_0`
  - \[electrons/m\] Initial MP size.
* - `Nx_regen`, `Ny_regen`, `Nvx_regen`, `Nvy_regen`, `Nvz_regen`
  - Number of grid points in each dimension of the uniform 5D phase-space grid used for
    regeneration.
* - `regen_hist_cut`
  - Phase-space density threshold below which no electrons are generated during regeneration.
* - `N_mp_soft_regen`
  - *Optional.* Soft regeneration threshold, checked at every bunch passage. Should be smaller
    than `N_mp_regen`. Soft regeneration is enabled only when this parameter and
    `N_mp_after_soft_regen` are set.
* - `N_mp_after_soft_regen`
  - *Optional.* Target number of MPs after a soft regeneration.
* - `N_mp_async_regen`, `N_mp_after_async_regen`
  - *Optional.* Threshold and target for asynchronous soft regeneration, which does not wait for
    the next bunch passage. Enabled only when both are set.
```

## Space charge

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `Dt_sc`
  - \[s\] Time step between updates of the electron space charge field map.
* - `Dh_sc`
  - \[m\] Grid size of the space charge particle-in-cell solver.
* - `t_sc_ON`
  - \[s\] Electron space charge forces are neglected for $t <$ `t_sc_ON`. Usually set to 0.
* - `sparse_solver`
  - *Optional, default `'PyKLU'`.* Sparse solver: `'PyKLU'` or `'scipy_slu'` (SciPy SuperLU).
    The old value `'klu'` is still accepted and selects `'PyKLU'`.
* - `PyPICmode`
  - *Optional, default `'FiniteDifferences_ShortleyWeller'`.* Field solver:
    `'FiniteDifferences_ShortleyWeller'`, `'ShortleyWeller_WithTelescopicGrids'`, or
    `'FiniteDifferences_Staircase'`.
* - `Dh_electric_energy`
  - *Optional, default `None`.* \[m\] Grid size used to compute the stored electric energy. If
    `None`, the electric energy is not computed.
* - `flag_em_tracking`
  - *Optional, default `False`.* Compute the forces from the electrons using electromagnetic
    potentials instead of the quasi-electrostatic approximation
    ([details](https://indico.cern.ch/event/840676/contributions/3532754/)).
* - `flag_reinterp_fields_at_substeps`
  - *Optional, default `False`.* Interpolate the beam and cloud field maps at every Boris
    sub-step.
```

## Multigrid

Used with `PyPICmode = 'ShortleyWeller_WithTelescopicGrids'`.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `f_telescope`
  - Magnification factor between successive grids. Must satisfy $0 < f < 1$.
* - `target_grid`
  - Target grid, given as a dictionary with the keys `'x_min_target'`, `'x_max_target'`,
    `'y_min_target'`, `'y_max_target'`, and `'Dh_target'`.
* - `N_nodes_discard`
  - Number of nodes at the edge of each internal grid that are discarded when interpolating the
    field (towards the coarser grid).
* - `N_min_Dh_main`
  - Minimum size of the first internal grid, in units of the main grid's $\Delta h$.
```

## Saving

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `filen_main_outp`
  - *Optional, default `'Pyecltest.mat'`.* Main output filename or path.
* - `save_only`
  - *Optional, default `None`.* Subset of output variables to save; all others are discarded.
* - `dec_fact_out`
  - Save the time-step output every `dec_fact_out` time steps.
* - `save_mat_every`
  - *Optional, default `1`.* Write the `.mat` output file every `save_mat_every` bunch passages.
* - `Dx_hist`
  - \[m\] Bin width of the horizontal histograms.
* - `r_center`
  - \[m\] Radius of the circle around the chamber centre (0, 0) used to compute the central
    electron density.
* - `Dt_En_hist`
  - \[s\] Time interval between saved energy spectra of the electrons hitting the chamber.
* - `Nbin_En_hist`
  - Number of bins of the energy histogram.
* - `En_hist_max`
  - \[eV\] Upper limit of the energy spectrum. Higher energies go into the last bin.
* - `flag_En_hist_seg`
  - *Optional, default `False`.* Save energy histograms for each chamber segment.
* - `flag_lifetime_hist`
  - *Optional, default `False`.* Save electron lifetime histograms. If `True`,
    `Dt_lifetime_hist`, `Nbin_lifetime_hist`, and `lifetime_hist_max` must also be set.
* - `Dt_lifetime_hist`
  - \[s\] Time interval between saved electron lifetime histograms.
* - `Nbin_lifetime_hist`
  - Number of bins of the electron lifetime histogram.
* - `lifetime_hist_max`
  - \[s\] Upper limit of the lifetime histogram. Longer lifetimes go into the last bin.
* - `flag_cos_angle_hist`
  - *Optional, default `True`.* Save a histogram of the cosine of the impact angle of electrons
    hitting the chamber.
* - `cos_angle_width`
  - *Optional, default `0.05`.* Bin width of the impact-angle histogram.
* - `flag_movie`
  - *Optional, default `0`.* Set to 1 to save the electron density distribution at every space
    charge evaluation.
* - `flag_sc_movie`
  - *Optional, default `0`.* Set to 1 to save the electron space charge field at every space
    charge evaluation.
* - `save_mp_state_time_file`
  - *Optional, default: no MP state is saved.* \[s\] Times at which MP positions and velocities
    are saved to file, given as a list or as a `.mat` file.
* - `flag_detailed_MP_info`
  - *Optional, default `0`.* Set to 1 to save the number of MPs at every time step.
* - `flag_hist_impact_seg`
  - *Optional, default `0`.* For polygonal chambers, set to 1 to save the number of electrons
    hitting each segment of the polygon.
* - `flag_verbose_file`, `flag_verbose_stdout`
  - *Optional, default `False`.* Report details of pathological impacts to an output file and to
    standard output, respectively.
* - `dec_fac_secbeam_prof`
  - *Optional, default `1`.* Decimation factor that reduces the number of points of the
    secondary beams' longitudinal profiles.
* - `el_density_probes`
  - *Optional, default: no probes.* Electron density probes, given as a list of Python
    dictionaries. Each probe needs a position (`x`, `y`) and a circle radius (`r_obs`).
* - `save_simulation_state_time_file`
  - *Optional, default: no state is saved.* \[s\] Times at which the full simulation state is
    saved to pickle files, given as a list or as a `.mat` file. See
    [Saving and restarting simulations](../tutorials/restart.md).
* - `checkpoint_folder`
  - *Optional, default: no checkpoints.* Folder in which checkpoints are saved.
* - `checkpoint_DT`
  - *Optional, default: no checkpoints.* \[s\] Simulated time between checkpoints.
* - `copy_main_outp_folder`
  - *Optional, default: the output is not copied.* Folder for backups of the output file. Backups
    are needed to restart a simulation from a checkpoint.
* - `copy_main_outp_DT`
  - *Optional, default: the output is copied only when a checkpoint is saved.* \[s\] Interval at
    which the output file is backed up.
* - `x_min_hist_det`, `x_max_hist_det`, `y_min_hist_det`, `y_max_hist_det`, `Dx_hist_det`
  - *Optional, default: no detailed histogram.* \[m\] Region of the chamber, and bin width, for
    an additional detailed horizontal histogram (`_det` stands for *detailed*).
* - `step_by_step_custom_observables`, `pass_by_pass_custom_observables`,
    `save_once_custom_observables`
  - *Optional, default: no custom output.* Custom observables saved at every time step, at every
    bunch passage, or once. See the examples for usage.
```

## Cloud particles

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `cloud_mass`
  - *Optional, default: electron mass.* \[kg\] Mass of the cloud particles.
* - `cloud_charge`
  - *Optional, default: electron charge.* \[C\] Charge of the cloud particles.
```

## Additional clouds

Simulations with several clouds are enabled with the following parameter. See
[Simulations with multiple clouds](multiple-clouds.md) for details.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `additional_clouds_file_list`
  - *Optional, default `[]`.* List of additional cloud file names.
```

## SEY curve extraction

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `extract_sey`
  - *Optional, default `True`.* Set to `False` to skip extracting the SEY curves.
```

## Emission energy distribution extraction

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `extract_ene_dist`
  - *Optional, default `False`.* Extract and save the energy distributions of secondary
    emission.
* - `ene_dist_test_E_impact_eV`
  - \[eV\] Impact energy used to extract the energy distributions.
* - `Nbin_extract_ene`
  - Number of bins of the extracted energy distributions.
* - `factor_ene_dist_max`
  - Sets the upper limit of the extracted energy distributions, which is
    `factor_ene_dist_max * ene_dist_test_E_impact_eV`.
```

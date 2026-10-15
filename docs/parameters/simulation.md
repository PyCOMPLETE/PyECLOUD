# Simulation parameters

This input file has to be named **simulation_parameters.input**.

**Other input filenames.** The following variables specify the names of the other three input
files which define the physical model of the simulation.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **machine_param_file**
  - Name of the machine parameter file.
* - **secondary_emission_parameters_file**
  - Name of the secondary emission configuration file.
* - **beam_parameters_file**
  - Name of the beam parameter file (in case of a multiple beam simulation this is the master
    beam).
```

**Secondary beams** The code can simulate the EC buildup in the presence of more than one
circulating beam. For this purpose a list of secondary beam files (one for each beam) has to be
provided. In the presence of secondary beams, the master beam determines the length of the
simulation, the energy used for the calculation of the bending field, and the bunch spacing
used for regeneration and saving purposes.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **secondary_beam_file_list**
  - (optional – default = \[\]) List of names of secondary beam parameter files.
```

**Log and progress files**

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **logfile_path**
  - Path of a text file that reports some synthetic information on the ongoing simulation (Passage
    number, number of MPs and number of electrons in the chamber).
* - **progress_path**
  - Path of a text file that reports the simulation progress in percent.
```

**Time sampling** The simulated time interval is defined by the length of the beam profile
(number of bunch passages) specified in the beam description.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **Dt**
  - \[s\] Simulation time step. This input is not considered if beam profile is imported from file.
* - **t_end**
  - \[s\] Extra time interval after the specified beam profile.

    Only coasting beam present in this interval. The value of this variable can also be negative in
    order to stop the simulation before the end of the specified profile.
```

**Negligible beam linear density**

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **lam_th**
  - \[beam part/m\] When the linear density of the beam is below this value, primary electrons are
    not generated and beam forces on the electrons are neglected.
```

**MP management settings**

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **N_mp_max**
  - Size of the arrays allocated to store the MP sizes, positions and velocities. MP regeneration
    settings must be set such that this number of MPs is never exceeded.
* - **N_mp_regen**
  - MP regeneration threshold. When the total number of MPs exceeds N_mp_regen a regeneration of
    the set of MPs is applied and a new Nref is chosen in order to get a target number of MPs. The
    check is performed at each t=n\*b_spac with n an integer.
* - **N_mp_regen_low**
  - Lower MP regeneration threshold. When the total number of MPs becomes smaller than
    N_mp_regen_low a regeneration of the set of MPs is applied and a new Nref is chosen in order to
    get a target number of MPs. The check is performed at each t=n\*b_spac with n an integer.
* - **t_ON_regen_low**
  - For t\< t_ON_regen_low the lower MP regeneration check is not performed.
* - **N_mp_after_regen**
  - Target number of MPs obtained after regeneration.
* - **fact_split**
  - MP split factor; when secondary emission occurs, if the size of the emitted charge is smaller
    than fact_split\*Nref, the impacting MP is simply rescaled according to the SEY. Otherwise new
    MPs are generated in order to keep the MP size as close as possible to Nref.
* - **fact_clean**
  - MP clean factor. MPs smaller than fact_clean\*Nref are deleted from the MP set. This test is
    performed at each t=n\*b_spac with n an integer.
* - **nel_mp_ref_0**
  - \[e-/m\] Initial MP size.
* - **Nx_regen, Ny_regen, Nvx_regen, Nvy_regen, Nvz_regen**
  - Number of mesh points in each dimension for the generation of the uniform grid for the 5-D
    phase space used for the MP regeneration.
* - **regen_hist_cut**
  - Defines a threshold value for the phase space density, below which no electrons are generated.
* - **N_mp_soft_regen**
  - (optional: by default soft regeneration is not applied, the presence of these inputs enables
    it)

    Soft regeneration threshold. The check is performed at each t=n\*b_spac with n an integer.
    N_mp_soft_regen should be smaller than N_mp_regen.
* - **N_mp_after_soft_regen**
  - (optional: by default soft regeneration is not applied, the presence of these inputs enables
    it)

    Target number of MPs obtained after a soft regeneration.
* - **N_mp_async_regen, N_mp_after_async_regen**
  - (optional: by default asynchronous soft regeneration is not applied, the presence of these
    inputs enables it)

    Parameters for asynchronous soft regeneration (the simulation does not wait for the next bunch
    passage to regenerate).
```

**Space charge parameters**

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **Dt_sc**
  - \[s\] Time step for the update of the electron space charge field map.
* - **Dh_sc**
  - \[m\] Grid size for the space charge PIC solver.
* - **t_sc_ON**
  - \[s\] Electron space charge forces are neglected for t\<t_sc_ON – normally t_sc_ON=0.
* - **sparse_solver**
  - (optional) Sparse solver: `PyKLU` (default) or `scipy_slu` (SciPy SuperLU).
* - **PyPICmode**
  - (optional) Options are ’FiniteDifferences_ShortleyWeller’ \[default\],
    ’ShortleyWeller_WithTelescopicGrids’, ’FiniteDifferences_Staircase’
* - **Dh_electric_energy**
  - \[m\] (optional – default = None) Grid size used to compute stored electric energy. If None is
    given, electric energy is not computed.
* - **flag_em_tracking**
  - (optional – default = False) Forces generated by electrons are computed using electromagnetic
    potentials instead of quasi-electrostatic approximation (for more information see
    <https://indico.cern.ch/event/840676/contributions/3532754/>).
* - **flag_reinterp_fields_at_substeps**
  - (optional – default = False) If true, field maps from beams and clouds are interpolated at each
    Boris substep.
```

**Multigrid parameters**

To be used with PyPICmode = ’ShortleyWeller_WithTelescopicGrids’

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **f_telescope**
  - Magnification factor between grids (it must be always $0<f<1$).
* - **target_grid**
  - Target grid parameters. It is a dictionary containing: ’x_min_target’, ’x_max_target’,
    ’y_min_target’, ’y_max_target’ and ’Dh_target’.
* - **N_nodes_discard**
  - Number of nodes at the edge of the internal grids which are discarded in the field
    interpolation (going to coarser grid).
* - **N_min_Dh_main**
  - Minimum size of the first internal grid in the $\Delta h$ of the main grid.
```

**Saving settings**

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **filen_main_outp**
  - Main output filename or path. (optional - default ’Pyecltest.mat’)
* - **save_only**
  - Define subset of output variables to be saved, other variables are discarded. (optional -
    default None)
* - **dec_fact_out**
  - Output is saved only for every nth time step.
* - **save_mat_every**
  - Output file .mat is saved only for every nth bunch passage. (optional - default 1)
* - **Dx_hist**
  - \[m\] Defines the binning for the saving of the horizontal histograms.
* - **r_center**
  - \[m\] Radius of a circle around the (0, 0) point (center of the chamber) for the calculation of
    the local central density.
* - **Dt_En_hist**
  - \[s\] Time step used to store the energy spectrum of the electrons impacting on the chamber.
* - **Nbin_En_hist**
  - Number of bins used for the energy histogram.
* - **En_hist_max**
  - \[eV\] Maximum energy in energy spectrum. Larger energies are attributed to the last bin of the
    histogram.
* - **flag_En_hist_seg**
  - (optional – default False) If enabled, energy histograms per segment are saved.
* - **flag_lifetime_hist**
  - If set to True enables lifetime histogram (optional – default=False).
```

If **flag_lifetime_hist** is set to True, the variables **Dt_lifetime_hist**,
**Nbin_lifetime_hist** and **lifetime_hist_max** need to be set.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **Dt_lifetime_hist**
  - \[s\] Time step used to store the histogram of the electrons lifetime.
* - **Nbin_lifetime_hist**
  - Number of bins used for the histogram of the electrons lifetime.
* - **lifetime_hist_max**
  - \[s\] Maximum lifetime in the histogram of the electrons lifetime. Larger lifetimes are
    attributed to the last bin of the histogram.
* - **flag_cos_angle_hist**
  - Save the cosine of the incident angle of impacting electrons. (optional – default True)
* - **cos_angle_width**
  - Width of the histogram bins. (optional – default 0.05)
* - **flag_movie**
  - (optional – default=0)

    (1 $\Rightarrow$ On, 0 $\Rightarrow$ Off) Saves files with the electron density
    distribution at each space charge evaluation.
* - **flag_sc_movie**
  - (optional– default=0)

    (1 $\Rightarrow$ On, 0 $\Rightarrow$ Off) Saves files with the electron space charge
    electric field distribution at each space charge evaluation.
* - **save_mp_state_time_file**
  - (optional – default: no MP state is saved)

    \[s\] List of instants at which MP positions and velocities are saved on file. The list can
    also be provided as a .mat file.
* - **flag_detailed_MP_info**
  - (optional – default=0)

    (1 $\Rightarrow$ On, 0 $\Rightarrow$ Off) Enables save of MP number at each time step.
* - **flag_hist_impact_seg**
  - (optional – default=0)

    (1 $\Rightarrow$ On, 0 $\Rightarrow$ Off) If the chamber is a polygon, enables saving the
    number of electrons which impacts onto each segment of polygon.
* - **flag_verbose_file, flag_verbose_stdout**
  - (optional – default=False)

    True/False toggles to save detailed information of pathological impacts in output file and
    stdout respectively.
* - **dec_fac_secbeam_prof**
  - (optional – default=1)

    Decimation factor to decrease the number of points of the longitudinal beam profile of
    secondary beams.
* - **el_density_probes**
  - (optional – default: no probes)

    (List of Python dictionaries) Probes for electron density evaluation. For each of them it is
    necessary to indicate the position (x,y) and radius of the circle (r_obs).
* - **save_simulation_state_time_file**
  - (optional – default: no simulation state is saved)

    \[s\] List of instants at which the full simulation state is saved on file (pickle files). The
    list can also be provided as a .mat file.
* - **checkpoint_folder**
  - (optional – default: no checkpoints are saved)

    Path to folder where checkpoints are saved.
* - **checkpoint_DT**
  - (optional – default: no checkpoints are saved)

    Time interval in seconds between checkpoints being saved.
* - **copy_main_outp_folder**
  - (optional – default: output is not copied)

    Path to folder for storing backups of the output. Backups are needed to restart a simulation
    from a checkpoint.
* - **copy_main_outp_DT**
  - (optional – default: output is only copied if a checkpoint is saved)

    Enables backing up the output file in regular intervals.
* - **x_min_hist_det, x_max_hist_det, y_min_hist_det, y_max_hist_det, Dx_hist_det**
  - (optional – default: no detailed histogram is saved)

    \[m\] Defines a region of the chamber where a detailed horizontal distribution histogram is
    saved (\_det stands for detailed).
* - **step_by_step_custom_observables, pass_by_pass_custom_observables,
    save_once_custom_observables**
  - (optional – default: no custom output)

    Custom observable definition. See example for detailed usage.
```

**Cloud particles**

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **cloud_mass**
  - (optional – default=electron mass) \[kg\] Mass of cloud particles.
* - **cloud_charge**
  - (optional – default=electron charge) \[C\] Charge of cloud particles.
```

**Additional clouds** Simulations with multiple clouds can be enabled with the following input
parameter. See the [multiple-clouds page](multiple-clouds.md) for a detailed description of this simulation mode.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **additional_clouds_file_list**
  - (optional – default=\[\]) List of names of additional cloud files.
```

**SEY** Extraction of the SEY curves can be disabled.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **extract_sey**
  - (optional – default=True) If set to False, it disables the extraction of the SEY curves.
```

**Extraction of electron emission energy distribution**

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **extract_ene_dist**
  - (optional – default=False) If set to True, secondary emission energy distributions are
    extracted and saved.
* - **ene_dist_test_E_impact_eV**
  - The impacting energy used in the energy distribution extraction. \[eV\]
* - **Nbin_extract_ene**
  - The number of bins used in the extracted energy distributions.
* - **factor_ene_dist_max**
  - Defines the maximum energy included in the extracted energy distribution. The maximum energy is
    factor_ene_dist_max \* ene_dist_test_E_impact_eV
```

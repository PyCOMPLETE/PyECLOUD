# Beam parameters

## Basic definitions

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `m0_part`
  - *Optional, default: proton mass.* \[kg\] Mass of the beam particles.
* - `q_part`
  - *Optional, default: proton charge.* \[C\] Charge of the beam particles.
* - `energy_eV`
  - \[eV\] Total energy of the beam particles.
* - `Dp_p`
  - *Optional, default `0`.* RMS momentum spread of the beam. Not used if the transverse beam
    size is given directly through `sigmax` and `sigmay`.
* - `nemittx`, `nemitty`
  - \[m\] Normalized transverse emittances of the beam. Not used if the transverse beam size is
    given directly through `sigmax` and `sigmay`.
* - `x_beam_pos`, `y_beam_pos`
  - *Optional, default `0`.* \[m\] Transverse position of the beam in the chamber.
* - `sigmax`, `sigmay`
  - *Optional, default `-1`.* \[m\] Transverse RMS beam size. If set to `-1`, the beam size is
    calculated from the energy, normalized emittance, dispersion, and momentum spread.
```

## Transverse electric field of the beam

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `beam_field_file`
  - Source of the transverse electric field map of the beam:

    - Name of a `.mat` file containing the field map.
    - `'computeFD'` (or `-1`): computed with the built-in finite-difference (FD) Poisson solver.
    - `'computeBE'`: computed with the Bassetti–Erskine formula plus image terms for the
      elliptical boundary. Only available for elliptical chambers.
    - `'compute_FDSW_multigrid'`: computed with the built-in FD Poisson solver on nested grids.
* - `save_beam_field_file_as`
  - *Optional, default: the field map is not saved.* Name of the file in which the field map is
    saved.
```

### Finite-difference solver (`'computeFD'`)

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `Dh_beam_field`
  - \[m\] Grid size of the FD solver and of the saved field map.
```

### Bassetti–Erskine (`'computeBE'`)

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `Nx`, `Ny`
  - Number of points of the field map in each plane.
* - `nimag`
  - Number of image terms.
```

### Multigrid solver (`'compute_FDSW_multigrid'`)

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `Dh_beam_field`
  - \[m\] Grid size of the FD solver and of the saved field map.
* - `f_telescope_beam`
  - Magnification factor between successive grids. Must satisfy $0 < f < 1$.
* - `target_grid_beam`
  - Target grid, given as a dictionary with the keys `'x_min_target'`, `'x_max_target'`,
    `'y_min_target'`, `'y_max_target'`, and `'Dh_target'`.
* - `N_nodes_discard_beam`
  - Number of nodes at the edge of each internal grid that are discarded when interpolating the
    field (towards the coarser grid).
* - `N_min_Dh_main_beam`
  - Minimum size of the first internal grid, in units of the main grid's $\Delta h$.
```

## Longitudinal beam profile

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `b_spac`
  - \[s\] Bunch spacing. Used to build the profile of a bunched beam. It is also the interval
    at which the MP set is cleaned and regenerated and at which data is saved.
* - `fact_beam`
  - Scaling factor applied to the beam profile.
* - `coast_dens`
  - *Optional, default `0`.* \[particles/m\] Coasting beam density added to the beam profile.
* - `flag_bunched_beam`
  - Set to 1 to generate the beam profile from a filling pattern, or to 0 to load it from a
    file.
```

### Profile from a filling pattern (`flag_bunched_beam = 1`)

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `sigmaz`
  - \[m\] RMS bunch length.
* - `t_offs`
  - \[s\] Time offset of the longitudinal profile. With `t_offs = 0`, only half of the first bunch
    is simulated.
* - `filling_pattern_file`
  - Name of a file with the intensities of the bunches, with zeros for empty slots. The pattern
    can also be given directly as a Python list, e.g. `filling_pattern_file = 4*(72*[1.]+8*[0.])`.
```

### Profile from a file (`flag_bunched_beam = 0`)

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `beam_long_prof_file`
  - Name of a `.mat` file containing the longitudinal beam profile. This file also sets the
    simulation time step.
```

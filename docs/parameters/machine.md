# Machine parameters

## Chamber profile

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `chamb_type`
  - *Optional, default `'ellip'`.* Shape of the chamber cross-section:

    - `'ellip'`: elliptical.
    - `'rect'`: rectangular.
    - `'polyg'` or `'polyg_cython'`: polygonal.
```

### Elliptical and rectangular chambers

Required when `chamb_type` is `'ellip'` or `'rect'`:

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `x_aper`, `y_aper`
  - \[m\] Horizontal and vertical semi-apertures of the chamber cross-section.
```

### Polygonal chambers

Required when `chamb_type` is `'polyg'`:

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `filename_chm`
  - Name of the file containing the horizontal and vertical coordinates of the chamber vertices.
    The vertices must be listed in counter-clockwise order; otherwise an error is raised.
* - `filename_chm_photoem`
  - Name of the file containing the chamber vertices used by the `'per_segment'` photoemission
    mode. The file must also contain the cumulative distribution function of the emission
    probability for each segment. These vertices must lie on the edges of the main chamber.
* - `flag_assume_convex`
  - *Optional, default `True`.* Assume the polygon is convex, which allows a faster
    impact-detection algorithm. An error is raised if the polygon is not convex; set to `False` to
    use the general (slower) algorithm for non-convex chambers.
```

## Tracking and magnetic field

Choose the tracking algorithm according to the magnetic field in the simulated element.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `track_method`
  - *Optional, default `'StrongBdip'`.* Tracking algorithm:

    - `'Boris'`: magnetic field map loaded from a file.
    - `'StrongBdip'`: uniform vertical magnetic field.
    - `'StrongBgen'`: generic transverse magnetic field map.
    - `'BorisMultipole'`: arbitrary multipoles, with the Boris tracker.
```

### `'Boris'`

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `B_map_file`
  - Name of a `.mat` file containing the magnetic field map (x, y, Bx(x, y), By(x, y),
    Bz(x, y)). For a quadrupole with unit gradient, set
    `B_map_file = 'analytic_quadrupole_unit_grad'`; no file is needed.
* - `fact_Bmap`
  - *Optional, default `1.`* Scaling factor applied to the field map.
* - `B0x`, `B0y`, `B0z`
  - *Optional, default `0`.* \[T\] Uniform field components added to the map.
```

### `'StrongBdip'`

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `B`
  - \[T\] Vertical magnetic field. If `B = -1`, the field is calculated from the beam energy and
    the total length of the bending magnets, `bm_totlen`.
* - `bm_totlen`
  - \[m\] Total length of the dipoles in the accelerator. Required when `B = -1`.
```

### `'StrongBgen'`

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `B_map_file`
  - Name of a `.mat` file containing the transverse field map. For a quadrupole with unit
    gradient, set `B_map_file = 'analytic_quadrupole_unit_grad'`; no file is needed.
* - `fact_Bmap`
  - *Optional, default `1.`* Scaling factor applied to the field map.
* - `B0x`, `B0y`
  - *Optional, default `0`.* \[T\] Uniform field components added to the map.
* - `B_zero_thrhld`
  - \[T\] Field magnitude below which the local magnetic field is treated as zero.
```

### `'BorisMultipole'`

`N_sub_steps` and `B_multip` are required; `B_skew` is optional.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `N_sub_steps`
  - Number of tracking sub-steps per time step.
* - `B_multip`
  - Normal multipole coefficients. Higher-order, skew, and combined-function magnets are
    supported. PyECLOUD uses the convention

    $B_y + iB_x = \sum_{n=0}^\infty \frac{1}{n!}\left( b_n + ib'_n \right)\left( x + iy \right)^{n}$

    where `B_multip` gives $b_n =\frac{\partial^n B_y}{\partial x^n}$. For example, a dipole with
    B = 8.3 T is `B_multip = [8.3]`, and a sextupole with 100 T/m² is
    `B_multip = [0., 0., 100.]`.
* - `B_skew`
  - *Optional, default `None`.* Skew multipole coefficients,
    $b'_n =\frac{\partial^n B_x}{\partial x^n}$. For example, a skew quadrupole with a gradient
    of 12 T/m is `B_multip = [0., 0.]` with `B_skew = [0., 12.]`.
* - `B0x`, `B0y`, `B0z`
  - *Optional, default `0`.* \[T\] Uniform field components added to the multipole field.
```

## Optics

These parameters can also be given in the beam parameter file(s), separately for each beam.
They are not needed if the transverse beam size is set directly in the beam parameter files.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `betafx`, `betafy`
  - \[m\] Horizontal and vertical beta functions at the simulated section.
* - `Dx`, `Dy`
  - *Optional, default `0`.* \[m\] Horizontal and vertical dispersion at the simulated section.
```

## Residual gas ionization

Primary electrons from residual gas ionization are generated only if these parameters are set.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `gas_ion_flag`
  - *Optional, default `0`.* Set to 1 to generate primary electrons by residual gas ionization.
* - `P_nTorr`
  - \[nTorr\] Pressure in the vacuum chamber.
* - `sigma_ion_MBarn`
  - \[MBarn\] Ionization cross-section of the residual gas.
* - `Temp_K`
  - \[K\] Temperature of the residual gas.
* - `unif_frac`
  - Fraction of primary electrons distributed uniformly in the chamber. The rest follow the
    transverse distribution of the beam.
* - `E_init_ion`
  - \[eV\] Initial energy of the generated electrons.
* - `t_ion`
  - \[s\] Electrons are generated by residual gas ionization only for $t <$ `t_ion`. The default
    (`-1`) generates none, so set `t_ion` (e.g. `t_ion = 100.`) when `gas_ion_flag = 1`. A
    warning is issued if `gas_ion_flag = 1` and `t_ion` ≤ 0.
```

## Photoemission

Primary electrons from photoemission are generated only if these parameters are set.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `photoem_flag`
  - *Optional, default `0`.* Photoemission mode:

    - `0`: off.
    - `1`: traditional model (see below).
    - `2` or `'from_file'`: angular distribution read from a file (see below).
    - `3` or `'per_segment'`: emission defined per chamber segment, through
      `filename_chm_photoem`.
* - `k_pe_st`
  - \[1/m\] Number of photoelectrons generated per beam particle and per unit length.
* - `energy_distribution`
  - *Optional, default `'gaussian'`.* Energy distribution of the photoelectrons, with parameters
    $\mu$ (`e_pe_max`) and $\sigma$ (`e_pe_sigma`):

    - `'gaussian'`: $p(E) = \frac{1}{\sigma\sqrt{2\pi}}\exp\left(-\frac{(E-\mu)^2}{2\sigma^2}\right)$
    - `'lognormal'`: $p(E) = \frac{1}{E\sigma\sqrt{2\pi}}\exp\left(-\frac{(\log E-\mu)^2}{2\sigma^2}\right)$
    - `'rect'`: $p(E)$ uniform for $\mu-\sigma/2 < E < \mu+\sigma/2$
    - `'mono'`: $E = \mu$
    - `'lorentz'`: $p(E) = \frac{1}{\pi}\frac{\sigma}{\sigma^2+(E-\mu)^2}$

    The Gaussian and Lorentz distributions are cut at zero so that no photoelectron has a
    negative energy.
* - `e_pe_sigma`
  - \[eV\] $\sigma$ of the energy distribution.
* - `e_pe_max`
  - \[eV\] $\mu$ of the energy distribution.
* - `out_radius`
  - \[m\] Radius of a circle enclosing the chamber. Photons are traced from inside the chamber
    towards points on this circle, and photoelectrons are generated where these paths cross the
    chamber wall. Must be larger than the chamber.
* - `phem_resc_fac`
  - *Optional, default `0.9999`.* Factor applied to the position of each generated
    photoelectron, so that it starts just inside the chamber wall.
* - `photoelectron_angle_distribution`
  - Must be set when `photoem_flag` is not 0. Emission angle distribution: `'cosine_3D'` or
    `'normal_emission'`. See `secondary_angle_distribution` in the
    [secondary emission parameters](secondary-emission.md).
```

### Traditional model (`photoem_flag = 1`)

Part of the photoelectrons is generated where unreflected photons hit the beam screen, and the
rest, from reflected photons, elsewhere on the chamber wall.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `inv_CDF_refl_photoem_file`
  - Name of a `.mat` file with the inverse cumulative distribution function (CDF) of the angular
    distribution of reflected photons. Set to `'unif_no_file'` for a uniform angular
    distribution; no file is then needed.
* - `refl_frac`
  - Fraction of photoelectrons generated by reflected photons.
* - `x0_refl`, `y0_refl`
  - \[m\] Impact point of the unreflected photons. `x0_refl` can also be set to `'left'` or
    `'right'`, in which case the point is computed automatically.
* - `alimit`
  - \[rad\] RMS width of the Gaussian angular distribution of photoelectrons generated by
    unreflected photons.
```

### Distribution from file (`photoem_flag = 2`)

The positions of all photoelectrons are drawn from a distribution read from a file.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `inv_CDF_all_photoem_file`
  - Like `inv_CDF_refl_photoem_file`, but for all photons. Set to `'unif_no_file'` for uniform
    generation over the whole chamber surface. Example scripts that create a file in the correct
    format are in `examples/photoemission_angular_distribution/` in the repository.
```

## Uniform initial distribution

The simulation starts with a given number of electrons distributed uniformly in the chamber.
Enabled only if these parameters are set.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `init_unif_flag`
  - *Optional, default `0`.* Set to 1 to start with a uniform electron distribution.
* - `Nel_init_unif`
  - \[electrons/m\] Initial number of electrons.
* - `E_init_unif`
  - *Optional, default `0`.* \[eV\] Initial energy of the electrons.
* - `x_max_init_unif`, `x_min_init_unif`, `y_max_init_unif`, `y_min_init_unif`
  - *Optional, default: the chamber semi-apertures, ±`x_aper` and ±`y_aper`.* \[m\] Edges of the
    region in which the electrons are generated.
```

## Uniform initial density

The simulation starts with a uniform electron volume density in a region of the chamber.
Enabled only if these parameters are set.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `init_unif_edens_flag`
  - *Optional, default `0`.* Set to 1 to start with a uniform electron density.
* - `init_unif_edens`
  - \[m$^{-3}$\] Initial electron volume density.
* - `E_init_unif_edens`
  - *Optional, default `0`.* \[eV\] Initial energy of the electrons.
* - `x_max_init_unif_edens`, `x_min_init_unif_edens`, `y_max_init_unif_edens`,
    `y_min_init_unif_edens`
  - *Optional, default: the chamber semi-apertures, ±`x_aper` and ±`y_aper`.* \[m\] Edges of the
    region in which the electrons are generated.
```

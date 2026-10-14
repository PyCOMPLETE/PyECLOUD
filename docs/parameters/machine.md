# Machine Parameters

**Chamber profile description**

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **chamb_type**
  - (optional – default = ‘ellip’)\
    Possible settings:

    1.  chamb_type = ‘ellip’ for elliptical chamber profile.
    2.  chamb_type = ‘rect’ for rectangular chamber profile.
    3.  chamb_type = ‘polyg’ or ‘polyg_cython’ for polygonal chamber profile.
```

When chamb_type = ‘ellip’ or ’rect’ the following input variables must be provided:

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **x_aper, y_aper**
  - Horizontal and vertical semi-apertures of the transverse chamber profile.
```

When chamb_type = ‘polyg’ the following input variable must be provided:

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **filename_chm**
  - Name of file containing horizontal and vertical vertexes of the chamber profile.
* - **filename_chm_photoem**
  - Name of file containing horizontal and vertical vertexes of the chamber profile that is used
    for the photoemission mode ’per_segment’. The chamber file must also contain the cumulative
    distribution function of the emission probability for each segment. The vertexes of this
    chamber must lie on the edges of the main chamber.
* - **flag_assume_convex**
  - \[optional\] Default: True
* - **flag_counter_clockwise_chamb**
  - \[optional\] Default: True. Only needed for the photoemission model ’per_segment’. It specifies
    the order in which the vertexes are defined.
```

**Tracking and magnetic field** (the tracking algorithm has to be chosen according to the
magnetic field conditions).

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **track_method**
  - (optional – default = ‘StrongBdip’)

    Possible settings:

    1.  track_method = ‘Boris’ for a magnetic field loaded from a file.
    2.  track_method = ‘StrongBdip’ for vertical uniform magnetic field.
    3.  track_method = ‘StrongBgen’ for a generic transverse magnetic field map.
    4.  track_method = ‘BorisMultipole’ for arbitrary multipoles with Boris electron tracker.
```

When track_method = ‘Boris’ the following input variables must be provided:

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **B_map_file**
  - Name of a .mat file containing the magnetic field map (x, y, Bx(x,y), By(x,y), Bz(x,y)). The
    case of a quadrupole magnet with unit gradient is embedded in the code and can be used setting:
    B_map_file = ‘analytic_quadrupole_unit_grad’.
* - **fact_Bmap**
  - (optional – default = 1.) Scaling factor applied to the magnetic field map.
* - **B0x, B0y, B0z**
  - (optional – default = 0) \[T\] Uniform magnetic fields added to the map.
```

When track_method = ‘StrongBdip’ the following input variables must be provided:

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **B**
  - \[T\] Value of the vertical uniform magnetic field. If B=-1 the magnetic field is calculated
    from the total length of the bending magnets (bm_totlen) and from the beam energy.
* - **bm_totlen**
  - \[m\] Total length of dipoles inside the accelerator (it must be provided when B=-1).
```

When track_method = ‘StrongBgen’ the following input variables must be provided:

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **B_map_file**
  - Name of a .mat file containing the transverse field map. The case of a quadrupole magnet with
    unit gradient is embedded in the code and can be used setting: B_map_file =
    ‘analytic_quadrupole_unit_grad’.
* - **fact_Bmap**
  - (optional – default = 1.) Scaling factor applied to the magnetic field map.
* - **B0x, B0y, B0z**
  - (optional – default = 0) \[T\] Uniform magnetic fields added to the map.
* - **B_zero_thrld**
  - \[T\] Value below which the local magnetic field is approximated with 0.
```

When track_method = ‘BorisMultipole’, N_sub_steps and B_multip must be provided, while B_skew
is optional:

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **N_sub_steps**
  - Number of tracking sub-steps per time step.
* - **B_multip**
  - Magnetic momenta. Higher order multipoles and skew magnets are supported, as well as compound
    magnets. The formalism in PyECLOUD follows this formula:
* - —
  - $B_y + iB_x = \sum_{n=0}^\infty \frac{1}{n!}\left( b_n + ib'_n \right)\left( x + iy \right)^{n}$
* - —
  - B_multip specifies $b_n =\frac{\partial^n B_y}{\partial x^n}$. For example: to simulate a
    dipole with B=8.3 T set B_multip = \[8.3\]; to simulate a sextupole with 100 T/m$^2$ set
    B_multip = \[0., 0., 100.\].
* - **B_skew**
  - The skew parameters can be included through B_skew. For a skew quadrupole with a gradient of 12
    T/m, set B_multip to \[0.,0.\] and B_skew to \[0.,12.\]. B_skew specifies
    $b'_n =\frac{\partial^n B_x}{\partial x^n}$ and defaults to *None*.
* - **B0x, B0y, B0z**
  - (optional – default = 0) \[T\] Uniform magnetic fields added to the map.
```

**Optics parameters** (can be provided also in the beam parameter file(s) independently for the
different beams) — not needed if transverse beam size is directly defined in the beam parameter
files.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **betafx, betafy**
  - \[m\] Horizontal and vertical beta functions at the simulation section.
* - **Dx, Dy**
  - (optional – default = 0) \[m\] Horizontal and vertical dispersion functions at the simulation
    section.
```

**Residual gas ionization parameters** (if the following input parameters are omitted primary
electron generation by residual gas ionization is not enabled).

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **gas_ion_flag**
  - (optional – default=0) (1 $\Rightarrow$ On, 0 $\Rightarrow$ Off) Enables primary electron
    generation by residual gas ionization.
* - **P_nTorr**
  - \[nTorr\] Pressure in the vacuum chamber.
* - **sigma_ion_MBarn**
  - \[MBarn\] Ionization cross section of the residual gas.
* - **Temp_K \[K\]**
  - Temperature of the residual gas.
* - **unif_frac**
  - Fraction of primary electrons uniformly distributed inside the chamber. The remaining part is
    generated according to the transverse distribution of the beam.
* - **E_init_ion**
  - \[eV\] Initial energy of the generated electrons.
* - **t_ion**
  - ???
```

**Photoemission parameters** (if generation by photoemission is not desired, the following
parameters can be omitted).

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **photoem_flag**
  - (optional – default=0) (0 $\Rightarrow$ Off, 1 $\Rightarrow$ traditional, 2 or ’from_file’
    $\Rightarrow$ Angular distribution from file, 3 or ’per_segment’ $\Rightarrow$ Part of
    chamber definition) Enables primary electron generation by photoemission.
* - **k_pe_st**
  - \[m$^{-1}$\] Number of photoelectrons to be generated per proton and per unit length.
* - **energy_distribution**
  - Can be one of:

    - ’gaussian’:
      $p(E) = \frac{1}{\sigma\sqrt{2\pi}}\exp\left(\frac{(E-\mu)^2}{2\sigma^2}\right)$
    - ’log-normal’: $p(E) = \frac{1}{E\sigma\sqrt{2\pi}}e^{\frac{(\log(E)-\mu)^2}{2\sigma^2}}$
    - ’rect’: $p(E)$ uniform for $\mu-\sigma < E < \mu+\sigma$
    - ’mono’: $E = \mu$
    - ’lorentz’: $p(E) = \frac{1}{\pi}\frac{\sigma}{\sigma^2+(E-\mu)^2}$

    The Lorentz and normal distributions are cut at 0 to prevent negative energies of photoemitted
    electrons.
* - **e_pe_sigma**
  - \[eV\] $\sigma$ from above.
* - **e_pe_max**
  - \[eV\] $\mu$ from above.
* - **out_radius**
  - \[m\] Radius of a circle external to the chamber.
* - **phem_resc_fact**
  - Rescaling factor on the position vector of the generated photoelectrons.
* - **photoelectron_angle_distribution**
  - See secondary_angle_distribution in the secondary emission parameters.
```

**Photoemission parameters for photoem_flag = 1**

**A part of photoelectrons is generated in the corner of the beam screen on direct impact of
photons (non-reflected). The reflected photons generated elsewhere.**

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **inv_CDF_refl_photoem_file**
  - Name of a .mat file providing the inverse of the Cumulative Distribution Function (CDF) for the
    angular distribution of the reflected photons. If inv_CDF_refl_photoem_file = ‘unif_no_file’
    the reflected photons have uniform angular distribution and no file needs to be provided.
* - **refl_frac**
  - Fraction of photoelectrons generated by reflected photons.
* - **x0_refl, y0_refl**
  - \[m\] Impact poiny of primary (not reflected) photons. x0_refl can be set to ‘left’ or ‘right’,
    in that case the point is computed automatically.
* - **alimit**
  - \[rad\] Extent (1$\sigma$) of the Gaussian angular distribution of photoelectrons generated
    by non-reflected photons.
```

**Photoemission parameters for photoem_flag = 2 or ’from_file’**

**The coordinates of all generated photoelectrons is specified from a file**

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **inv_CDF_all_photoem_file**
  - See inv_CDF_refl_photoem_file, but for all photons. If it is set to ’unif_no_file’, the
    photoelectron generation is uniform over the whole surface. Example scripts to generate a file
    with the correct format can be found in the folder
    `examples/photoemission_angular_distribution/` in the repository.
```

**Uniform initial distribution** Simulation starts with electrons uniformly distributed in the
chamber (if the following input parameters are omitted this feature is not enabled).

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **init_unif_flag**
  - (optional – default=0) (1 $\Rightarrow$ On, 0 $\Rightarrow$ Off) Enables uniform electron
    distribution at the beginning of the simulation.
* - **Nel_init_unif**
  - Initial number of electrons.
* - **E_init_unif**
  - \[eV\] Initial energy of the generated electrons.
* - **x_max_init_unif, x_min_init_unif, y_max_init_unif, y_min_init_unif**
  - \[m\] Edges of the region where the electrons are generated.
```

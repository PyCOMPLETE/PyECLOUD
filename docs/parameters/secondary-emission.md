# Secondary emission parameters

## Choice of the model

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `switch_model`
  - *Optional, default `0`.* Secondary emission model:

    - `0` or `'ECLOUD'`: the model used in the ECLOUD code.
    - `1` or `'ACC_LOW'`: a more accurate treatment of low-energy impacts.
    - `'ECLOUD_nunif'`: the ECLOUD model, with SEY parameters given per segment in the chamber
      file.
    - `'ECLOUD_nunif_charging'`: the ECLOUD model with charging effects, with SEY and charging
      parameters given per segment in the chamber file.
    - `'cos_low_ene'`: cosine dependence at low energy.
    - `'flat_low_ene'`: flat dependence at low energy.
    - `'furman_pivi'`: the Furman–Pivi model.
    - `'from_file'`: interpolation of a user-supplied SEY curve (see `sey_file`).
    - `'perfect_absorber'`: every impacting electron is absorbed; no secondaries are emitted.
```

## ECLOUD model

These parameters are described in G. Iadarola's
[doctoral thesis](https://repository.cern/records/pb0hw-d9008).

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `Emax`
  - \[eV\] Impact energy at which the SEY is maximal.
* - `del_max`
  - Maximum SEY.
* - `s_param`
  - *Optional; only for `switch_model = 'ECLOUD'`.* Shape parameter of the true secondary SEY
    curve.
* - `R0`
  - Weight of the reflected (elastic) electron component.
* - `E0`
  - *Optional, default set by the model.* \[eV\] Shape parameter of the reflected SEY curve.
* - `E_th`
  - \[eV\] Maximum energy of true secondary electrons.
* - `sigmafit`
  - $\sigma$ of the log-normal energy distribution of true secondaries.
* - `mufit`
  - $\mu$ of the log-normal energy distribution of true secondaries.
```

## Furman–Pivi model

The symbols refer to the Furman–Pivi paper (M. A. Furman and M. T. F. Pivi, Phys. Rev. ST Accel.
Beams 5, 124404, 2002). See also
[CERN-ACC-NOTE-2019-0029](https://cds.cern.ch/record/2683285).

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `E_th`
  - \[eV\] Maximum energy of true secondary electrons.
* - `furman_pivi_surface`
  - Python dictionary with all Furman–Pivi model parameters. Its keys are listed below.
```

### General keys

```{list-table}
:header-rows: 1
:widths: 35 65

* - Key
  - Description
* - `'use_modified_sigmaE'`
  - bool. Modify $\sigma_e$ in the same way as POSINST does.
* - `'use_ECLOUD_theta0_dependence'`
  - bool. Use the same incidence-angle dependence of the SEY as the ECLOUD model.
* - `'use_ECLOUD_energy'`
  - bool. Use the same secondary energy distribution as the ECLOUD model. Elastically reflected
    electrons then also use the ECLOUD emission angle.
* - `'conserve_energy'`
  - bool. Ensure that the total energy of the electrons emitted in one event does not exceed the
    energy of the impacting electron.
* - `'exclude_rediffused'`
  - bool. Set $\delta_r = 0$.
* - `'choice'`
  - `'poisson'` or `'binomial'`: probability distribution used for $P^\prime_{n,ts}$.
* - `'M_cut'`
  - int. Maximum number of macroparticles emitted in a single event.
* - `'p_n'`
  - NumPy array of `M_cut` floats: $p_n$.
* - `'eps_n'`
  - NumPy array of `M_cut` floats: $\epsilon_n$.
```

### Backscattered electrons

```{list-table}
:header-rows: 1
:widths: 35 65

* - Key
  - Description
* - `'p1EInf'`
  - float. $P_{1,e}(\infty)$.
* - `'p1Ehat'`
  - float. $\hat{P}_{1,e}$.
* - `'eEHat'`
  - float. $\hat{E}_e$.
* - `'w'`
  - float. $W$.
* - `'p'`
  - float. $p$.
* - `'e1'`
  - float. $e_1$.
* - `'e2'`
  - float. $e_2$.
* - `'sigmaE'`
  - float. $\sigma_e$.
* - `'theta_e_max'`
  - *Optional; by default the angle is not limited.* float, \[rad\]. Upper limit applied to the
    incidence angle when computing the angular dependence for backscattered electrons. Must
    satisfy $0 < \theta_{e,\max} < \pi/2$.
```

### Rediffused electrons

```{list-table}
:header-rows: 1
:widths: 35 65

* - Key
  - Description
* - `'p1RInf'`
  - float. $P_{1,r}(\infty)$.
* - `'eR'`
  - float. $E_r$.
* - `'r'`
  - float. $r$.
* - `'q'`
  - float. $q$.
* - `'r1'`
  - float. $r_1$.
* - `'r2'`
  - float. $r_2$.
* - `'theta_r_max'`
  - *Optional; by default the angle is not limited.* float, \[rad\]. Upper limit applied to the
    incidence angle when computing the angular dependence for rediffused electrons. Must
    satisfy $0 < \theta_{r,\max} < \pi/2$.
```

### True secondaries

```{list-table}
:header-rows: 1
:widths: 35 65

* - Key
  - Description
* - `'deltaTSHat'`
  - float. $\hat{\delta}_{ts}$.
* - `'eHat0'`
  - float. $\hat{E}_{ts}$.
* - `'s'`
  - float. $s$.
* - `'t1'`
  - float. $t_1$.
* - `'t2'`
  - float. $t_2$.
* - `'t3'`
  - float. $t_3$.
* - `'t4'`
  - float. $t_4$.
```

## Other parameters

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `sey_file`
  - Required for `switch_model = 'from_file'`. Path to a file giving the reflected and true
    secondary yields as a function of energy. Examples are in `examples/sey_files/` in the
    repository.
* - `flag_costheta_delta_scale`
  - *Optional, default `True`.* Scale the true secondary SEY curve with the impact angle.
* - `flag_costheta_Emax_shift`
  - *Optional, default `True`.* Shift the true secondary SEY curve in energy with the impact
    angle.
* - `switch_no_increase_energy`
  - *Optional, default `0`.*

    - `0`: off.
    - `1`: no longer used.
    - `2`: ensure that true secondaries do not have more energy than the impacting electron.
* - `thresh_low_energy`
  - \[eV\] Maximum impact energy for which the `switch_no_increase_energy` check is applied.
* - `scrub_en_th`
  - \[eV\] Minimum impact energy of scrubbing electrons, used to estimate the scrubbing current.
* - `secondary_angle_distribution`
  - Required. Emission angle distribution of new electrons, where $\theta$ is the angle between
    the surface normal and the initial velocity:

    - `'cosine_3D'`: $dn/d\theta = \cos\theta\sin\theta$, which includes the surface element
      $\sin\theta$ in spherical coordinates
      ([reference](http://www.sciencedirect.com/science/article/pii/S0042207X02001732)).
    - `'normal_emission'`: $\theta = 0$.

    The legacy `'cosine_2D'` distribution ($dn/d\theta = \cos\theta$) is no longer supported
    from PyECLOUD v8.9.0 and raises an error. To reproduce the behaviour of PyECLOUD 6.6.0 or
    earlier, use a version older than v8.9.0.
```

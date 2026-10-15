# Secondary emission parameters

**Choice of the model**

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **switch_model**
  - (optional – default=0)

    Different secondary emission models:

    1.  switch_model = 0 or ‘ECLOUD’ for the model used in ECLOUD.
    2.  switch_model = ‘ACC_LOW’ for a more accurate treatment of low energy impacts.
    3.  switch_model = ‘ECLOUD_nunif’ for the model used in ECLOUD, with sey info in the chamber
        shape file.
    4.  switch_model = ‘ECLOUD_nunif_charging’ for the model used in ECLOUD and charging effects,
        with sey and charging info in the chamber shape file.
    5.  switch_model = ‘cos_low_ene’ for cosine low energy dependence.
    6.  switch_model = ‘flat_low_ene’ for flat low energy dependence.
    7.  switch_model = ‘furman_pivi’ for the Furman-Pivi model of secondary electron emission.
    8.  switch_model = ‘from_file’ for interpolation of a user-specified curve.
```

**Secondary Electron Yield** These are the parameters for the ‘ECLOUD’ model, as described in
G. Iadarola’s thesis.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **Emax**
  - Energy corresponding to the maximum SEY.
* - **s_param**
  - Shape parameter of the true SEY curve.
* - **del_max**
  - Maximum of the SEY curve.
* - **R0**
  - Weight of the reflected electron component.
* - **E0**
  - Shape parameter of the reflected SEY curve.
* - **E_th**
  - Maximum energy for true secondary electrons.
* - **sigmafit**
  - Sigma parameter of the lognormal distribution.
* - **mufit**
  - Mu parameter of the lognormal distribution.
```

**Secondary Electron Yield** These are the parameters for the ‘furman_pivi’ model.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **E_th**
  - Maximum energy for true secondary electrons.
* - **furman_pivi_surface**
  - A Python dict containing all parameters needed for the Furman-Pivi model. The parameters are
    listed and described below.
* - General parameters. To be specified in furman_pivi_surface.
  - 1.  ’use_modified_sigmaE’: Boolean, modifies the $\sigma_e$ parameter from the FP paper in
        the same way as done in POSINST.
    2.  ’use_ECLOUD_theta0_dependence’: Boolean, use same incident angle dependence of SEY as in
        ECLOUD.
    3.  ’use_ECLOUD_energy’: Boolean, use same secondary electron energy distribution as in ECLOUD.
        Also makes elastic electrons use the ECLOUD emission angle.
    4.  ’conserve_energy’: Boolean, ensures the sum of emitted electron energies in a single event
        is not greater than the impacting electron energy.
    5.  ’exclude_rediffused’: Boolean, sets $\delta_r = 0$.
    6.  ’choice’: ’poisson’ or ’binomial’, choice of probability distribution for
        $P^\prime_{n,ts}$
    7.  ’M_cut’: int, the maximum number of macroparticles that can be emitted in a single event.
    8.  ’p_n’: numpy array containing M_cut floats. Corresponds to $p_n$ in FP paper.
    9.  ’eps_n’: numpy array containing M_cut floats. Corresponds to $\epsilon_n$ in FP paper.
* - Parameters for backscattered electrons. To be specified in furman_pivi_surface.
  - 1.  ’p1EInf’: float, $P_{1,e}(\infty)$.
    2.  ’p1Ehat’: float, $\hat{P}_{1,e}$.
    3.  ’eEHat’: float, $\hat{E}_e$.
    4.  ’w’: float, $W$.
    5.  ’p’: float, $p$.
    6.  ’e1’: float, $e_1$.
    7.  ’e2’: float, $e_2$.
    8.  ’sigmaE’: float, $\sigma_e$.
    9.  ’theta_e_max’: float, optional, \[rad\] value to which the incidence angle is limited
        when computing the angular dependence for backscattered electrons. Must satisfy
        $0 < \theta_{e,\max} < \pi/2$. If omitted, the angle is not limited.
* - Parameters for rediffused electrons. To be specified in furman_pivi_surface.
  - 1.  ’p1RInf’: float, $P_{1,r}(\infty)$.
    2.  ’eR’: float, $E_r$.
    3.  ’r’: float, $r$
    4.  ’q’: float, $q$.
    5.  ’r1’: float, $r_1$.
    6.  ’r2’: float, $r_2$.
    7.  ’theta_r_max’: float, optional, \[rad\] value to which the incidence angle is limited
        when computing the angular dependence for rediffused electrons. Must satisfy
        $0 < \theta_{r,\max} < \pi/2$. If omitted, the angle is not limited.
* - Parameters for true secondaries. To be specified in furman_pivi_surface.
  - 1.  ’deltaTSHat’: float, $\hat{\delta}_{ts}$.
    2.  ’eHat0’: float, $\hat{E}_{ts}$.
    3.  ’s’: float, $s$.
    4.  ’t1’: float, $t_1$
    5.  ’t2’: float, $t_2$
    6.  ’t3’: float, $t_3$
    7.  ’t4’: float, $t_4$
```

**Other parameters**

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **sey_file**
  - Needed for secondary emission model ‘from_file’. Path to a file that specifies the reflected
    and true secondary emission yields for different energies. An example can be found in the
    subfolder sey_files in the PyECLOUD code.
* - **flag_costheta_delta_scale**
  - If enabled the SEY curve associated to true secondary electrons is scaled as a function of the
    angle of the impacting electron. Default is TRUE.
* - **flag_costheta_Emax_shift**
  - If enabled the SEY curve associated to true secondary electrons is shifted (change of energy
    scale) as a function of the angle of the impacting electron. Default is TRUE.
* - **switch_no_increase_energy**
  - (Default=0)

    1.  switch_no_increase_energy = 0 : Off.
    2.  switch_no_increase_energy = 1 : No longer in use.
    3.  switch_no_increase_energy = 2 : Ensure that the true secondaries do not have more energy
        than the corresponding impacting electrons.
* - **thresh_low_energy**
  - Maximum energy for which the energy check is performed.
* - **scrub_en_th**
  - Minimum energy of scrubbing electrons (for scrubbing current estimations).
* - **secondary_angle_distribution**
  - (no default – must be set) Can be ’cosine_3D’ or ’normal_emission’. For new electrons,
    $\theta$ describes the angle between the surface normal and their initial velocity vector.
    For ’cosine_3D’, d$n$/d$\theta=\cos\theta\sin\theta$, which takes into account the surface
    element $\sin\theta$ in spherical coordinates. More info
    [here](http://www.sciencedirect.com/science/article/pii/S0042207X02001732). For
    ’normal_emission’, $\theta = 0$.

    The legacy ’cosine_2D’ distribution (d$n$/d$\theta=\cos\theta$) is no longer supported as of
    PyECLOUD v8.9.0 and raises an error. To reproduce the behaviour of PyECLOUD 6.6.0 or earlier,
    use a version older than v8.9.0.
```

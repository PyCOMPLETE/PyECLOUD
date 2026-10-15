# Buildup simulation outputs

The main output is a MATLAB `.mat` file, named `Pyecltest.mat` by default.
`filen_main_outp` selects another name. Multiple-cloud simulations write an
output file for each cloud. This page covers the buildup output; single-bunch
beam monitors use HDF5 files as shown in the
[single-bunch tutorial](../tutorials/single-bunch.md).

## Read a result

```python
from PyECLOUD.myloadmat_to_obj import myloadmat_to_obj

data = myloadmat_to_obj("Pyecltest.mat")
print(sorted(name for name in vars(data) if not name.startswith("__")))
print(data.t.shape, data.Nel_timep.shape)
```

The loader exposes output variables as attributes and removes singleton dimensions
by default (`squeeze=True`). Inspect shapes before indexing
histograms, especially for runs with only one sample or bunch passage.
Use `myloadmat_to_obj("Pyecltest.mat", squeeze=False)` to preserve singleton axes.
The file can also be read with SciPy's `scipy.io.loadmat`.
The [buildup tutorial](../tutorials/buildup.md#inspect-your-result) plots the
simulation time series.
See [Plotting buildup results](../tutorials/buildup-outputs.md) for the full
set of example figures and histogram projections.

## Time series

`N_steps` denotes the number of saved time samples, `N_pass` the number of
recorded bunch passages, `N_bins` the number of horizontal histogram bins,
and `N_en` the number of energy bins. The time-series arrays below have shape
`(N_steps,)` after squeezing, provided there is more than one sample.

Use the saved `t` array as the time coordinate. The saver uses a sampling
interval based on `Dt_ref * dec_fact_out`; saving and tracking need not occur
at the same cadence. Entries for emitted and impacting electrons accumulate
over the interval since the previous saved sample.

PyECLOUD represents a transverse section. Electron population sums are
normalized per unit longitudinal length, not counts integrated over a whole
magnet. Weighted electron energies likewise carry that longitudinal
normalization. Central and probe densities divide by the transverse sampling
area to give volume density.

| Variable | Meaning |
| --- | --- |
| `t` | Saved simulation times [s]; one value per recorded sample |
| `lam_t_array` | Beam line density [particles/m] at those times |
| `Nel_timep` | Electron population in the chamber [electrons/m] |
| `Nel_emit_time` | Electrons emitted during the saved interval, per metre |
| `Nel_imp_time` | Electrons impacting during the saved interval, per metre |
| `En_emit_eV_time` | Weighted emitted-electron energy during the interval [eV/m] |
| `En_imp_eV_time` | Weighted impacting-electron energy during the interval [eV/m] |
| `En_kin_eV_time` | Total weighted electron kinetic energy currently in the chamber [eV/m] |
| `cen_density` | Electron volume density in a transverse disk of radius `r_center` around the origin [m⁻³] |

These time-series arrays share the saved sample index. With
`flag_detailed_MP_info`, `N_mp_time` records the total number of macroparticles
at each saved time and `nel_mp_ref_time` records their reference weight.
Other enabled models can add fields, such as electric-field energy or
cross-ionization observables.

## Bunch passages and histograms

`t_hist` indexes bunch-passage records and `xg_hist` gives horizontal bin
coordinates. Horizontal histograms are organized by passage and horizontal
bin. Energy spectra use their own recorded times, `t_En_hist`, and energy-bin
coordinates, `En_g_hist`; do not assume their number of records equals the
number of bunch passages.

The passage arrays `t_hist`, `N_mp_pass`, `N_mp_ref_pass`, `N_mp_impact_pass`,
and `N_mp_corrected_pass` have shape `(N_pass,)`. The horizontal coordinates
`xg_hist` have shape `(N_bins,)`. The four horizontal histograms—`nel_hist`,
`nel_impact_hist_tot`, `nel_impact_hist_scrub`, and `energ_eV_impact_hist`—have
shape `(N_pass, N_bins)`, with one row per passage and one column per bin.
Singleton axes are removed by the loader's default squeezing.

| Variable | Meaning |
| --- | --- |
| `t_hist` | Times associated with bunch-passage records [s] |
| `xg_hist` | Horizontal histogram coordinates [m] |
| `N_mp_pass` | Number of macroparticles in the chamber at each passage |
| `N_mp_ref_pass` | Reference macroparticle weight at each passage |
| `N_mp_impact_pass` | Number of macroparticle impacts during each passage |
| `N_mp_corrected_pass` | Number of macroparticle tracking corrections during each passage |
| `nel_hist` | Horizontal histogram of the electron population |
| `nel_impact_hist_tot` | Horizontal histogram of electron impacts |
| `nel_impact_hist_scrub` | Horizontal impact histogram restricted to energies above `scrub_en_th` |
| `energ_eV_impact_hist` | Horizontal histogram of weighted energy deposited at the walls |
| `En_hist` | Impacting-electron counts per recorded interval and energy bin |
| `all_Ekin_hist` | Kinetic-energy histogram of all electrons in the chamber at each recorded energy-sampling time |
| `t_En_hist` | Times associated with the energy spectra [s] |
| `En_g_hist` | Energy-bin centres [eV] |
| `lifetime_hist` | Electron lifetime counts per recorded interval and lifetime bin, with `flag_lifetime_hist` |
| `lifetime_g_hist` | Lifetime-bin centres [s] |
| `t_lifetime_hist` | Saved lifetime-histogram times [s], controlled by `Dt_lifetime_hist` |

Lifetime histograms have their own sampling interval, `Dt_lifetime_hist`;
use `t_lifetime_hist` to locate their records in time.

The final bin of each energy and lifetime histogram collects entries outside
the selected energy or lifetime range. Counts
are binned accumulations; converting them into a density per unit energy or
time requires the corresponding bin width and accumulation interval.

## Other and optional fields

Some fields are absent when their feature is disabled; others may contain an
empty array or a sentinel value. Inspect the saved data before plotting them.

| Variable | Meaning |
| --- | --- |
| `dec_fact_out` | Time-series output sampling factor: for a fixed step `Dt_ref`, output is saved every `dec_fact_out` steps |
| `b_spac` | Bunch spacing used internally [s] |
| `t_sc_video` | Recorded space-charge update times [s] |
| `U_sc_eV` | Legacy field described as unused in the original manual |
| `chamber_area` | Transverse chamber area [m²] |
| `nel_hist_impact_seg` | Electron impacts on chamber segments, with `flag_hist_impact_seg` |
| `nel_hist_emit_seg` | Emitted-electron histogram on chamber segments |
| `energ_eV_impact_seg` | Total weighted energy deposited by impacting electrons on each chamber segment [eV/m], with `flag_hist_impact_seg` |
| `t_sec_beams` | Times of saved secondary-beam information [s] |
| `sec_beam_profiles` | Secondary-beam line-density profiles [particles/m] |
| `el_dens_at_probes` | Volume densities at locations configured by `el_density_probes` [m⁻³], indexed by probe and saved sample |
| `x_el_dens_probes`, `y_el_dens_probes` | Probe centre coordinates [m] |
| `r_el_dens_probes` | Probe radii [m] |
| `nel_hist_det` | Detailed electron-position histogram, with `flag_hist_det` |
| `xg_hist_det` | Coordinates of the detailed histogram bins [m], with `flag_hist_det` |
| `cos_angle_hist` | Histogram of the cosine of electron impact angles |
| `xg_hist_cos_angle` | Coordinates of the impact-angle histogram bins |

This reference adapts the output section of
{download}`the original manual <../legacy_tex/reference/reference.pdf>` by
Giovanni Iadarola, Eleonora Belli, Philipp Dijkstal, Lotta Mether, Annalisa Romano,
Giovanni Rumolo, and Eric Wulff. Model-specific fields are not exhaustively
listed here.

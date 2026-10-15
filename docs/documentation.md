# Working on the documentation

The website uses Sphinx, MyST Markdown, and the Read the Docs theme. Sources
live in `docs/`. Edit Markdown (`.md`) pages and add
new pages to a `toctree` so they appear in the navigation. The logo is
`docs/pyecl_logo_darkbg.png`; colours and layout adjustments live in
`docs/_static/custom.css`.

Parameter tables use MyST's fenced `list-table` directive so their Markdown
cells can contain paragraphs, lists, and equations without long table rows.

## Build and preview locally

From the repository root, create a separate environment in `docs/.venv-docs`
using Python 3.11 or newer (Read the Docs builds with Python 3.12):

```bash
python -m venv docs/.venv-docs
source docs/.venv-docs/bin/activate
python -m pip install -r docs/requirements.txt
python -m sphinx -n -W --keep-going -b html docs docs/_build/html
python -m http.server 8000 --directory docs/_build/html
```

Open <http://localhost:8000>. On Windows, activate the environment with
`docs\.venv-docs\Scripts\activate` instead. The environment is ignored by Git
and excluded from the Sphinx sources in `docs/conf.py`. Rebuild after editing a page, then refresh
the browser. Stop the server with Ctrl+C.

The build treats warnings as errors and checks internal references. It does
not install PyECLOUD, compile native extensions, or run simulations. The
displayed release is read from `PyECLOUD/_version.py` without importing the
simulation packages.

Check external links separately, since remote sites can be temporarily
unavailable:

```bash
python -m sphinx -b linkcheck docs docs/_build/linkcheck
```

## GitHub and Read the Docs

Run the `Documentation` GitHub Actions workflow manually from the **Actions**
tab. Download the `documentation-html` artifact from a successful run,
extract it, and serve the extracted folder with `python -m http.server 8000`.
Actions may need to be enabled in a fork's **Actions** tab first.

For a hosted preview, connect your fork to a separate Read the Docs project and
select the documentation branch in its version settings. The repository's
`.readthedocs.yaml` installs the same documentation requirements and enables
strict Sphinx builds. Enable pull request builds in the project's settings to
preview pull requests opened against your fork. This requires configuring the
hosting integration; committing the configuration alone does not publish a site.

## Updating content

The buildup parameter pages were migrated from
`docs/legacy_tex/reference/src/reference.tex`. The single-bunch page includes
`PyPARIS_sim_class/Simulation_parameters_doc.md` directly, so that guide has one
editable source. Preserve original authorship when migrating further material.

Check parameter names, units, defaults, and examples against the implementation
when changing their descriptions.

### Audit and validation notes

A clean documentation build checks markup and references; it does not verify
simulation behaviour. The tutorials still need a complete numerical validation.

The buildup parameter pages have been checked against
`PyECLOUD/default_input_parameters.py` and the code that reads each parameter:

- Every documented parameter name, and every `furman_pivi_surface` key, exists
  in the code. Misspelled names inherited from the manual were corrected
  (`secondary_beams_file_list`, `B_zero_thrhld`, `phem_resc_fac`,
  `theta_r_max`), and `flag_counter_clockwise_chamb`, which does not exist,
  was removed.
- Defaults were checked for parameters marked *Optional*, and accepted option
  values were checked for `switch_model`, `secondary_angle_distribution`,
  `energy_distribution`, `beam_field_file`, `sparse_solver`, and `track_method`.
- Behaviour described in the descriptions was checked where it was changed or
  added, for example `t_ion`, `fact_clean`, the uniform initial density, and
  the photoemission energy distributions.

Not yet checked: the physical meaning and units of the SEY model parameters,
most units elsewhere, and descriptions carried over unchanged from the manual.
The input parameters `stopfile`, `filename_init_MP_state`,
`flag_continuous_emission`, and `cloud_name` are not yet documented.

## Wiki migration record

The initial migration used the public PyECLOUD wiki at revision
`5d30c514549031a5d5baef3eba706a42fa27feed`, reviewed on 7 October 2026.
New pages link to their original sources so authorship and revision history
remain traceable.

The manual was converted to Markdown during the initial migration. Its
converter is not part of the documentation build; edit the Markdown sources
directly.

The restart guide was checked against `PyECLOUD/buildup_simulation.py`,
`PyECLOUD/pyecloud_saver.py`, and `PyPARIS_sim_class/Save_Load_Status.py`.
Output field names were checked against `PyECLOUD/pyecloud_saver.py`, including
`xg_hist_det`, which was duplicated as `nel_hist_det` in the original manual.
The output reference also preserves the manual's array dimensions, sampling
definitions, and enabling flags. The manual describes lifetime histograms as
being saved per bunch passage; the current implementation uses
`Dt_lifetime_hist` and saves `t_lifetime_hist`, so the reference documents that
timing instead. Energy normalization and the transverse sampling region for
`cen_density` follow the implementation rather than the manual's older wording.
The parallel execution guide retains the wiki's MPI environment advice;
the old Python and MPI source-compilation recipes remain historical.

The buildup tutorial plots the output of the reader's simulation. The older
plotting script in `docs/legacy_tex/example/` reads a stored regression reference file.
The plotting guide presents all 20 legacy output examples in five grouped SVG
figures. They are generated from the saved LHC dipole reference data using
`docs/_static/plot_buildup_outputs.py` and stored in `docs/_static/buildup/`.
Sphinx includes these figures directly, so building the website does not require
running PyECLOUD or regenerating plots. To regenerate them from the repository
root with PyECLOUD and Matplotlib installed:

```bash
python docs/_static/plot_buildup_outputs.py testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns/Pyecltest_angle3D_ref.mat --output-dir docs/_static/buildup
```

| Original material | Destination |
| --- | --- |
| Home: simulation modes, libraries, citation | Homepage, installation, and citation pages |
| Buildup tutorial | Updated buildup tutorial and plotting example |
| PyECLOUD–PyHEADTAIL tutorial | Single-bunch tutorial and parallel execution |
| Reference manual | Input parameter pages, output reference, and downloadable PDF |
| Physical models and algorithms | Physics references |
| Python 3 setup: MPI advice | Parallel execution guide |
| Miniforge installation | Link to the official installation guide |
| Checkpointing presentation link | Restart guide checked against the current code |

Python 2, Ubuntu 14.04, the Python 3 transition announcement, and the obsolete
index remain historical. The old macOS recipe has not been promoted to current
instructions without platform validation. The separate Furman–Pivi and fast
beam-ion wiki pages were placeholders; the useful Furman–Pivi references were
retained from the physics page. The old PDF compilation page is superseded by
the Sphinx build instructions; its LaTeX sources remain in `docs/legacy_tex/reference/src/`.

# Working on the documentation

The website uses Sphinx, MyST Markdown, and the Read the Docs theme. Sources
live in `docs/`. Edit Markdown (`.md`) pages and add
new pages to a `toctree` so they appear in the navigation. The logos are in
`docs/logos/`; the site uses `pyecloud-logo-stacked-transparent-for-dark.svg`
on the dark blue panel. Colours, taken from the logos, and layout adjustments
live in `docs/_static/custom.css`.

Parameter tables use MyST's fenced `list-table` directive so their Markdown
cells can contain paragraphs, lists, and equations without long table rows.

## Build and preview locally

From the repository root, create a separate environment in `docs/.venv-docs`
using Python 3.11 or newer (the publishing workflow uses Python 3.12):

```bash
python -m venv docs/.venv-docs
source docs/.venv-docs/bin/activate
python -m pip install -r docs/requirements.txt
python -m sphinx -n -W --keep-going -b html docs docs/_build/html
python -m http.server 8000 --directory docs/_build/html
```

Open [http://localhost:8000](http://localhost:8000). On Windows, activate the environment with
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

## Publishing

The site is published with GitHub Pages by the `Documentation` workflow
(`.github/workflows/docs.yml`):

- Every push to `main`, including merged pull requests, builds the site and
  deploys it.
- Every pull request against `main` builds the site without deploying it, so
  warnings and broken references are caught before merging.
- The workflow can also be started manually from the **Actions** tab.

The build uses the same strict settings as the local build above. Only
`PyCOMPLETE/PyECLOUD` deploys; forks run the build check only.

To preview a build, download the `documentation-html` artifact from the run
summary, extract it, and serve the extracted folder with
`python -m http.server 8000`. In a fork, Actions may need to be enabled in the
fork's **Actions** tab first.

## Updating content

The buildup parameter pages were migrated from
`docs/legacy_tex/reference/src/reference.tex`. The single-bunch page includes
`PyPARIS_sim_class/Simulation_parameters_doc.md` directly, so that guide has one
editable source. Parameter descriptions should match the code; check names, defaults and units in `PyECLOUD/default_input_parameters.py` when you change them.

### Buildup output figures

The figures in [Plotting buildup results](tutorials/buildup-outputs.md) are
generated from the saved LHC dipole reference data using
`docs/_static/plot_buildup_outputs.py` and stored in `docs/_static/buildup/`.
Sphinx includes these figures directly, so building the website does not require
running PyECLOUD or regenerating plots. To regenerate them from the repository
root with PyECLOUD and Matplotlib installed:

```bash
python docs/_static/plot_buildup_outputs.py testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns/Pyecltest_angle3D_ref.mat --output-dir docs/_static/buildup
```

---

How the documentation was migrated from the original wiki and reference manual,
and how far it has been checked against the code, is recorded in the
[migration and audit notes](migration-notes.md).

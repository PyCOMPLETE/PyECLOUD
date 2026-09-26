# PyECLOUD

PyECLOUD is a 2D macro-particle code for the simulation of electron cloud effects in particle accelerators.

## Installation

Python 3.11 or newer and working C and Fortran compilers are required. In a conda
environment, these can be installed with `conda install -c conda-forge c-compiler
fortran-compiler`. Activate the environment before building.

```sh
python -m pip install pyecloud
```

Pip automatically installs the Poisson solvers from **pypic-poisson**. Only source
distributions are published for these two packages, so pip builds their native
extensions locally. Their Python import names remain `PyECLOUD` and `PyPIC`.

### Development installation

With sibling PyPIC and PyECLOUD checkouts, run from their parent directory:

```sh
python -m pip install -e ./PyPIC
python -m pip install -e ./PyECLOUD
```

Or, after installing PyPIC, run `python -m pip install -e .` inside this checkout.
Pip installs the Python build dependencies and compiles seven Fortran extensions
with F2PY/Meson and two Cython/C extensions. There is no separate `make` or
`cythonize` step. Cython, Meson, and Ninja are build dependencies; NumPy, SciPy,
matplotlib, and pypic-poisson are runtime dependencies. The unrelated `pypic`
distribution on PyPI is not a dependency. If you previously installed this
PyPIC checkout under the old distribution name `PyPIC`, uninstall it before
reinstalling the renamed package to avoid overlapping installed files.

Python edits take effect immediately with an editable installation. Rerun the
installation command after changing native sources. Use `python -m pip install .`
for a regular installation. `make`, `setup_pyecloud`, and `cythonize` remain
convenience wrappers around the same editable pip installation.

Optional integrations can be installed with `.[pyheadtail]` (PyHEADTAIL and h5py)
or `.[hdf5]` (h5py). PyKLU is optional; the SciPy sparse solver is available with
the core dependencies.

## Running simulations

The Python namespace is unchanged:

```python
from PyECLOUD.buildup_simulation import BuildupSimulation

sim = BuildupSimulation(pyecl_input_folder="/path/to/input_folder")
sim.run()
```

The input folder contains `simulation_parameters.input`, machine and secondary
emission parameters, and beam files. Existing configuration/data paths retain
their original meaning; choose your working directory and output paths as before.

The launch scripts now live under `examples/`:

```sh
python examples/000_run_simulation.py /path/to/input_folder
python examples/001_reload_state_and_run.py /path/to/input_folder /path/to/simulation_state_0.pkl
```

## Repository layout and validation

- `PyECLOUD/`: Python modules, Cython/C sources, and `fortran/` sources.
- `tests/`: automated installation and numerical smoke tests.
- `testing/`: existing simulation regression cases and their reference data.
- `examples/`: simulation launch scripts.
- `other/`, `doc/`, and `dev/`: studies, documentation, and maintenance scripts.

```sh
python -m pip install -e '.[tests]'
python -m pytest
python -m pip install build
python -m build
```

The tests exercise the installed extensions and a short build-up simulation.
`python -m build` creates a source distribution and builds a wheel from it.
The large historical regression datasets stay in the repository and are not
included in the installed package. Version metadata lives in `PyECLOUD/_version.py`;
simulation logs include Git provenance when running from a checkout and work
without Git in a wheel installation.

## Publishing a source release

Publish the required `pypic-poisson` release first. Set the PyECLOUD version in
`PyECLOUD/_version.py`, then install the release tools and check the source archive:

```sh
python -m pip install build twine
python release.py --build-only
```

After testing the archive, commit and push the release changes, then run:

```sh
python release.py
```

The script requires a clean checkout and an unused `v<version>` tag. It builds
and checks one source archive, uploads only that `.tar.gz` to PyPI using your
Twine credentials (for example, configured in `~/.pypirc`), then creates and
pushes the version tag to `origin`. No wheels are uploaded. `--build-only`
retains the archive in `dist/` without uploading or tagging. If the upload
succeeds but tagging or pushing fails, finish those Git operations manually;
PyPI does not allow re-uploading the same release file.

More information about installation and usage can be found in the [Wiki](https://github.com/PyCOMPLETE/PyECLOUD/wiki).

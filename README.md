# PyECLOUD

PyECLOUD is a 2D macro-particle code for the simulation of electron cloud effects in particle accelerators.

## Installation

Create a Conda environment with Python 3.13 and the compilers and numerical
libraries, then install the released PyHEADTAIL and PyECLOUD packages:

```sh
conda create -n ecloud -c conda-forge python=3.13 pip
conda activate ecloud
conda install -c conda-forge c-compiler cxx-compiler fortran-compiler
conda install -c conda-forge numpy scipy matplotlib
pip install --upgrade pyheadtail
pip install --upgrade pyecloud
```

Keep the environment active when installing packages and running simulations.
See the [installation guide](docs/installation.md) for editable Git installs.

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
matplotlib, pypic-poisson, PyKLU, and h5py are runtime dependencies. The unrelated `pypic`
distribution on PyPI is not a dependency. If you previously installed this
PyPIC checkout under the old distribution name `PyPIC`, uninstall it before
reinstalling the renamed package to avoid overlapping installed files.

Python edits take effect immediately with an editable installation. Rerun the
installation command after changing native sources. Use `python -m pip install .`
for a regular installation. `make`, `setup_pyecloud`, and `cythonize` remain
convenience wrappers around the same editable pip installation.

PyHEADTAIL can also be installed through the `.[pyheadtail]` extra; h5py is a
core dependency. PyKLU is installed automatically as the default sparse solver
for its speed advantage over SciPy's SuperLU. PyECLOUD can also run without
PyKLU: select SciPy with `sparse_solver = 'scipy_slu'` in
`simulation_parameters.input`.

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

The launch scripts live under `examples/getting_started/`:

```sh
python examples/getting_started/000_run_simulation.py /path/to/input_folder
python examples/getting_started/001_reload_state_and_run.py /path/to/input_folder /path/to/simulation_state_0.pkl
```

## Repository layout and validation

- `PyECLOUD/`: Python modules, Cython/C sources, and `fortran/` sources.
- `tests/`: automated installation and numerical smoke tests.
- `testing/`: existing simulation regression cases and their reference data.
- `examples/`: launch scripts, example simulations, and studies.
- `docs/`: website sources and the original manual under `legacy_tex/`.
- `dev/`: maintenance scripts.

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

Publish the required `pypic-poisson` release first. Use a development environment
with the build dependencies installed. Set a new version in `PyECLOUD/_version.py`,
then build and check the source archive:

```sh
python -m pip install build twine
python -m build --sdist
python -m twine check --strict dist/*.tar.gz
```

After testing the archive, commit and push the release changes, then run from
the repository root:

```sh
./release.sh
```

The script reads the version from `PyECLOUD/_version.py`, builds the sdist with
`python -m build`, creates and pushes `v<version>` to `origin`, then uploads the
archive using Twine and your PyPI credentials (for example, configured in
`~/.pypirc`). The archive is built in a temporary directory so only this
release's source archive is uploaded, then cleaned up. No wheels are uploaded.
A failed build stops the script before anything is tagged. If uploading fails
after the tag is pushed, finish the upload manually.

## Documentation

The documentation site includes installation, simulation tutorials, input and
output references, and links to the physical models.
See [Working on the documentation](docs/documentation.md) for local builds,
previews, and publishing. Building the site does not require installing the
simulation packages or compiling their native extensions.

The original [Wiki](https://github.com/PyCOMPLETE/PyECLOUD/wiki) remains available
as a historical source. See the documentation's migration record for coverage.

<h1 align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/PyCOMPLETE/PyECLOUD/main/docs/logos/pyecloud-logo-inline-transparent-for-dark.svg">
    <img src="https://raw.githubusercontent.com/PyCOMPLETE/PyECLOUD/main/docs/logos/pyecloud-logo-inline-transparent-for-light.svg" alt="PyECLOUD" width="420">
  </picture>
</h1>

<p align="center">
  <a href="https://pypi.org/project/PyECLOUD/"><img src="https://img.shields.io/pypi/v/PyECLOUD?color=0A1330&label=PyPI" alt="PyPI version"></a>
  <img src="https://img.shields.io/badge/python-%E2%89%A53.11-0A1330" alt="Python 3.11 or newer">
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS-0A1330" alt="Platforms: Linux and macOS">
  <a href="https://github.com/PyCOMPLETE/PyECLOUD/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-Apache--2.0-0A1330" alt="License: Apache-2.0"></a>
  <a href="https://github.com/PyCOMPLETE/PyECLOUD/actions/workflows/docs.yml"><img src="https://github.com/PyCOMPLETE/PyECLOUD/actions/workflows/docs.yml/badge.svg?branch=main" alt="Documentation build status"></a>
</p>

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

PyECLOUD and its Poisson solver package, `pypic-poisson`, are published as
source distributions only, so pip compiles their native extensions during
installation; this is why the compilers are needed. Pip installs
`pypic-poisson` and the other dependencies automatically. The packages are
imported as `PyECLOUD` and `PyPIC`.

PyKLU is installed as the default sparse solver because it is faster than
SciPy's SuperLU. To use SciPy instead, set `sparse_solver = 'scipy_slu'` in
`simulation_parameters.input`.

PyHEADTAIL is only needed for beam-dynamics simulations. It can also be
installed together with the optional frequency-analysis package `nafflib`
through the `tracking` extra: `pip install "pyecloud[tracking]"`.

### Development installation

For development, install PyECLOUD in editable mode from a Git checkout. If you
are also working on PyPIC, clone both repositories side by side and run from
their parent directory:

```sh
python -m pip install -e ./PyPIC
python -m pip install -e ./PyECLOUD
```

Otherwise, run `python -m pip install -e .` inside the PyECLOUD checkout; pip
then installs `pypic-poisson` from PyPI. See the
[installation guide](docs/installation.md) for details.

Pip installs the build dependencies (Cython, Meson, and Ninja) and compiles
seven Fortran extensions with F2PY and two Cython/C extensions. Changes to Python
files take effect immediately. After changing Fortran, Cython, or C sources,
run the installation command again. `make` and `./setup_pyecloud` are
shortcuts for `python -m pip install -e .`.

Do not install the unrelated `pypic` package from PyPI; PyECLOUD's Poisson
solvers come from `pypic-poisson`.

## Running simulations

A buildup simulation reads its settings from an input folder:

```python
from PyECLOUD.buildup_simulation import BuildupSimulation

sim = BuildupSimulation(pyecl_input_folder="/path/to/input_folder")
sim.run()
```

The input folder contains `simulation_parameters.input`, the machine and
secondary emission parameter files, and the beam file. The output files are
written relative to the current working directory.

The same can be done with the scripts in `examples/getting_started/`, which
also show how to resume from a saved simulation state:

```sh
python examples/getting_started/000_run_simulation.py /path/to/input_folder
python examples/getting_started/001_reload_state_and_run.py /path/to/input_folder /path/to/simulation_state_0.pkl
```

## Repository layout and tests

- `PyECLOUD/`: Python modules, Cython/C sources, and `fortran/` sources.
- `PyPARIS/`, `PyPARIS_sim_class/`, `PyPARIS_CoupledBunch_sim_class/`: launchers
  and simulation classes for coupled PyECLOUD–PyHEADTAIL simulations.
- `tests/`: automated installation and numerical smoke tests.
- `testing/`: simulation regression cases and their reference data. They are
  kept in the repository and are not part of the installed package.
- `examples/`: launch scripts, example simulations, and studies.
- `docs/`: documentation sources, and the original reference manual under
  `legacy_tex/`.

To run the tests and build the package:

```sh
python -m pip install -e '.[tests]'
python -m pytest
python -m pip install build
python -m build
```

The tests check the compiled extensions and run a short buildup simulation.
`python -m build` creates a source distribution and builds a wheel from it.

The package version is set in `PyECLOUD/_version.py`. Each simulation log
records the PyECLOUD version and, when running from a Git checkout, the commit.

## Publishing a source release

Publish the required `pypic-poisson` release first. Set the new version in
`PyECLOUD/_version.py`, then build and check the source archive:

```sh
python -m pip install build twine
python -m build --sdist
python -m twine check --strict dist/*.tar.gz
```

After testing the archive, commit and push the release changes to `main`, then
run from the repository root:

```sh
./release.sh
```

The script reads the version from `PyECLOUD/_version.py`, checks that the
current commit is on `main` and that the sdist builds and passes `twine check`,
then creates and pushes the tag `v<version>`. A failed check stops the script
before anything is tagged. The tag is pushed to `origin`; set `REMOTE` to push
elsewhere, e.g. `REMOTE=upstream ./release.sh`.

Pushing the tag to `PyCOMPLETE/PyECLOUD` starts the **Publish to PyPI**
workflow (`.github/workflows/publish.yml`). It checks that the tagged commit is
on `main` and that the tag matches `_version.py`, builds and checks the sdist,
and uploads it to PyPI using trusted publishing, so no PyPI token is needed.
No wheels are uploaded. Follow the run in the repository's **Actions** tab.
Tags pushed to forks build the archive but do not publish.

## Documentation

The documentation site includes installation, simulation tutorials, input and
output references, and links to the physical models.
See [Working on the documentation](docs/documentation.md) for local builds,
previews, and publishing. Building the site does not require installing the
simulation packages or compiling their native extensions.

The original [Wiki](https://github.com/PyCOMPLETE/PyECLOUD/wiki) remains available
as a historical source. See the [migration and audit notes](docs/migration-notes.md)
for how its content was carried over.

# Electron cloud buildup

A buildup simulation tracks the electron cloud at a section of an accelerator.
The prescribed beam generates fields and primary electrons; the cloud does not
change the beam trajectory in this mode.

## Obtain the example inputs

Complete the [installation](../installation.md) first. The examples and reference
datasets are kept in Git; a PyPI installation does not supply the full example
tree. If you do not already have a checkout, clone the main PyECLOUD repository:

```bash
git clone https://github.com/PyCOMPLETE/PyECLOUD.git
cd PyECLOUD
```

If you already have a checkout, start from its root instead. You can use the
installed PyPI package with these inputs; an editable install is only needed
when developing the code. Use inputs from a revision compatible with your
installed release.

From the checkout root, copy the LHC dipole example into a new working directory:

```bash
PYECL_SOURCE="$PWD"
PYECL_CASE="$PYECL_SOURCE/testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns"
mkdir ../ecloud-buildup
cp "$PYECL_CASE"/*.input ../ecloud-buildup/
cp "$PYECL_CASE/beam.beam" "$PYECL_CASE/LHC_chm_ver.mat" ../ecloud-buildup/
cd ../ecloud-buildup
```

## Configure and run

The working directory contains:

| File | Purpose |
| --- | --- |
| `simulation_parameters.input` | Numerical settings, output settings, and names of other input files |
| `machine_parameters.input` | Chamber, magnetic field, and primary-electron generation |
| `secondary_emission_parameters.input` | Secondary emission model |
| `beam.beam` | Beam profile |
| `LHC_chm_ver.mat` | Chamber geometry used by the example |

See the [input reference](../parameters/index.md) before changing these settings.
This example uses a large particle population and can take substantial
computing time.

Save the following as `run_buildup.py` in the working directory:

```python
from PyECLOUD.buildup_simulation import BuildupSimulation

simulation = BuildupSimulation(pyecl_input_folder=".")
simulation.run()
```

Run it with the environment active:

```bash
python run_buildup.py
```

The repository also provides `examples/getting_started/000_run_simulation.py`, which accepts an
input-folder argument. Running from the example's working directory keeps its
relative input and output paths together.

## Inspect your result

The default main output is `Pyecltest.mat`. Save this as `plot_buildup.py` beside
your output file:

```python
import matplotlib.pyplot as plt
from PyECLOUD.myloadmat_to_obj import myloadmat_to_obj

data = myloadmat_to_obj("Pyecltest.mat")
fig, axes = plt.subplots(2, 1, sharex=True)
axes[0].plot(data.t * 1e9, data.lam_t_array)
axes[0].set_ylabel("Beam particles / m")
axes[1].plot(data.t * 1e9, data.Nel_timep)
axes[1].set_ylabel("Electrons / m")
axes[1].set_xlabel("Time [ns]")
fig.tight_layout()
fig.savefig("buildup.png", dpi=150)
```

```bash
python plot_buildup.py
```

The resulting `buildup.png` shows the beam profile and cloud population.

For more examples, see [Plotting buildup results](buildup-outputs.md).
# Simulations with multiple clouds

A simulation can include more than one cloud. Each cloud feels the space charge of the others,
and cross-ionization processes between clouds can be defined.

## Main input

Add these parameters to `simulation_parameters.input`:

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `additional_clouds_file_list`
  - *Optional, default `[]`.* List of additional cloud file names.
* - `cross_ion_definitions`
  - *Optional, default `None`.* Dictionary describing the cross-ionization processes (see this
    [presentation](https://indico.cern.ch/event/835473/contributions/3525107/) for details).
```

## Cloud files

Each cloud file contains a subset of the parameters of the main input files.

### Required parameters

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - `cloud_mass`, `cloud_charge`
  - See [Simulation parameters](simulation.md#cloud-particles).
* - `gas_ion_flag`, `photoem_flag`, `init_unif_flag`, `init_unif_edens_flag`
  - See [Machine parameters](machine.md).
* - `switch_model`
  - See [Secondary emission parameters](secondary-emission.md).
```

### Optional parameters

A cloud file can also contain the following parameters. Any parameter missing from the cloud
file takes its value from the main input files, or its default if it is not set there either.

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameters
  - Reference
* - MP management, saving, log and progress files
  - [Simulation parameters](simulation.md)
* - `N_sub_steps`
  - [Machine parameters: `'BorisMultipole'`](machine.md#borismultipole)
* - Residual gas ionization, photoemission, uniform initial distribution and density
  - [Machine parameters](machine.md)
* - Secondary emission
  - [Secondary emission parameters](secondary-emission.md)
```

Each cloud writes its own `.mat` output file, named after the main output file and the cloud
file. The contents are described in the [output reference](../reference/outputs.md).

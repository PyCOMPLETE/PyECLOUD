# Simulations with multiple clouds

The code provides the option of including more than one cloud in the simulation. Each cloud is
influenced by the other clouds through their space charge. Cross-ionization processes between
different clouds can be defined. These options can be enabled by providing the following
additional information.

In the **simulation_parameters.input** file:

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **additional_clouds_file_list**
  - (optional – default=\[\]) List of names of additional cloud files.
* - **cross_ion_definitions**
  - (optional – default={ }) Dictionary that describes the desired cross-ionization processes (more
    details can be found in <https://indico.cern.ch/event/835473/contributions/3525107/>).
```

Each cloud file should contain a subset of the input parameters specified in the main input
files, as listed below. The following parameters are required in the cloud input file.

**Mandatory cloud parameters**

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **cloud_mass**, **cloud_charge**
  - See the [Simulation parameters](simulation.md) section.
* - **gas_ion_flag**, **photoem_flag**, **init_unif_flag**, **init_unif_edens_flag**
  - See the [Machine parameters](machine.md) section.
* - **switch_model**
  - See the [Secondary emission parameters](secondary-emission.md) section.
```

In addition, a cloud file can contain the following optional parameters. Any optional parameter
that is not specified in the cloud input file takes its value from the main input, and gets the
default value if the parameter is not given in the main input files.

**Optional cloud parameters**

```{list-table}
:header-rows: 1
:widths: 35 65

* - Parameter
  - Description
* - **MP management settings**, **Saving settings**, **Log and progress files**
  - See the [Simulation parameters](simulation.md) section for the detailed list of parameters
    and their usage.
* - **N_sub_steps**
  - See under “Tracking and magnetic field” in the [Machine parameters](machine.md) section.
* - **Residual gas ionization parameters**, **Photoemission parameters**, **Uniform initial
    distribution**, **Uniform initial density**
  - See the [Machine parameters](machine.md) section for the detailed list of parameters and
    their usage.
* - **Secondary emission parameters**
  - As listed in the [Secondary emission parameters](secondary-emission.md) section.
```

Each cloud produces its own .mat output file, as described in the
[output reference](../reference/outputs.md).
The file is named according to the main output file name and the cloud file name.

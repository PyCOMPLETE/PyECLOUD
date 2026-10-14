# Saving and restarting simulations

PyECLOUD buildup simulations support explicit saved states and rolling
checkpoints. PyPARIS single-bunch simulations use a separate multijob status
and bunch-state mechanism.

## Save a buildup state at chosen times

In `simulation_parameters.input`, specify simulation times in seconds:

```python
save_simulation_state_time_file = [1e-6, 2e-6]
```

Choose times within the simulated interval. This enables state files such as
`simulation_state_0.pkl`. To load one, use the original input files and a
compatible software environment. From the simulation working directory:

```python
from PyECLOUD.buildup_simulation import BuildupSimulation

simulation = BuildupSimulation(pyecl_input_folder=".")
simulation.load_state("simulation_state_0.pkl")
simulation.run()
```

By default, `load_state` disables further explicit state saving and changes the
main output name to a restarted variant. The repository's
`examples/getting_started/001_reload_state_and_run.py` provides the same workflow with input-folder
and state-file arguments. These pickle files contain Python objects, so preserve
the environment used to produce them; they are not a portable interchange format.

## Rolling checkpoints for buildup runs

Set these inputs before starting the run:

```python
checkpoint_DT = 1e-6
checkpoint_folder = './checkpoints/'
```

`checkpoint_DT` is an interval in **simulated seconds**, not wall-clock seconds.
Use the trailing slash in the folder path because the current implementation
concatenates it with checkpoint filenames. Select an interval appropriate to
the simulated duration and checkpoint size.

Restart the original buildup command from the same working directory, keeping
the input files, main output, and checkpoint folder together. Initialization
automatically loads a checkpoint when the folder contains exactly one file. An
empty folder starts a fresh simulation; more than one file raises an error.
The saver normally removes the previous checkpoint after writing its replacement.

Checkpoint recovery restores the saved state and trims output recorded after
the checkpoint. Preserve a copy of your results before investigating a failed
recovery. Checkpointing does not submit another scheduler job by itself.

## Single-bunch multijob runs

`PyPARIS_sim_class` records progress in `simulation_status.sta` and saves the
bunch between job parts. Keep the status and bunch files together. A normal
subsequent launch continues after a completed part; an unfinished part is
rejected rather than silently skipped.

Configure `N_turns`, `N_turns_target`, `check_for_resubmit`, `submission_system`,
and, where needed, `resubmit_command` as described in the
[multijob parameter reference](../parameters/single-bunch.md#multijob-setup).
Keep `check_for_resubmit = False` for manual continuation. Do not remove the
status file as a way to resume: removing it starts the numbering from scratch.

Further background: [checkpointing presentation](https://indico.cern.ch/event/772318/).

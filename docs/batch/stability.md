# Stability simulations with continuation

For single-bunch PyPARIS simulations, use a DAGMan node that runs one part at
a time. The worker script downloads the inputs and last completed state, runs
one simulation part, and uploads its checkpoint. Read the
[storage and environment setup](index.md) first.

## Configure simulation parts

In `Simulation_parameters.py`, for example:

```python
N_turns = 100
N_turns_target = 1000
check_for_resubmit = True
submission_system = 'HTCondor'
resubmit_command = None
```

Choose each part to finish comfortably within the job time limit, including
input and output transfers. PyPARIS returns `177` when a completed part needs
continuation. A normal final exit is `0`; inspect the simulation logs to
distinguish reaching the target from an intentional early stop. Other codes
represent failures. See [saving and restarting](../tutorials/restart.md).

## Submit description

Save as `htcondor.sub`. Keep this file, `pyparis_executable.sh`, the POST
script, DAG, and `logs/` on AFS. Put the simulation parameters and
`pyecloud_config` in an existing EOS run directory. This example transfers
`LHC_chm_ver.mat` separately; replace its URL with your geometry file and make
sure the simulation configuration refers to the transferred filename.

```text
universe = vanilla
executable = pyparis_executable.sh
arguments = $(request_cpus)
output = logs/$(ClusterId).$(ProcId).out
error = logs/$(ClusterId).$(ProcId).err
log = logs/$(ClusterId).log
should_transfer_files = YES
transfer_input_files = root://eosproject-e.cern.ch//eos/project/e/ecloud-simulations/YOURDIRECTORY/stability/LHC_chm_ver.mat
when_to_transfer_output = ON_EXIT
transfer_output_files = ""
request_cpus = 8
request_memory = 16GB
request_disk = 20GB
+JobFlavour = "tomorrow"
queue 1
```

Adjust CPU, memory, and disk requests to your simulation.
Passing `request_cpus` to the script makes PyPARIS use the allocated number
of processes. Simulation data travels to EOS from
the worker, so no simulation filenames appear in `transfer_output_files`.

## Worker script

Download {download}`pyparis_executable.sh <../_static/htcondor/pyparis_executable.sh>`
and save it beside `htcondor.sub`. Set `ROOT_URL`, `SIM_PATH`, and the EOS
environment path. The script normalizes trailing slashes in `ROOT_URL` and
`SIM_PATH`. See
[Python environment](index.md#python-environment) for setup instructions.

This script follows the example in
`examples/HTCondor_templates/PyPARIS_singlebunch_checkpoint_job/`.

```{literalinclude} ../_static/htcondor/pyparis_executable.sh
:language: bash
```

### Checkpoint files

The same file set is downloaded on continuation and uploaded after each
completed part. `NN` is `present_simulation_part` formatted with at least two
digits (`00`, `01`, ...).

| File | Purpose |
| --- | --- |
| `simulation_status.sta` | Part number and completion state; uploaded last |
| `bunch_status_partNN.h5` | Saved bunch state |
| `bunch_evolution_NN.h5`, `slice_evolution_NN.h5` | Evolution monitors for the part |
| `sim_param.pkl` | Saved simulation parameters |
| `pyparislog.txt` | Simulation log |
| `envinfo.txt`, `stdout.txt`, `stderr.txt` | Environment information and worker-captured Python output |
| `multigrid_config_dip.pkl`, `multigrid_config_dip.txt` | Optional saved dipole multigrid configuration, copied when present |

`Simulation_parameters.py` and `pyecloud_config` are fetched for every part.
Keep these inputs and the software environment consistent throughout the run.
The run directory must exist before the first job. A successful directory
listing without a status file means a fresh run; a failed listing aborts the
job. Never remove the status file to recover a failed continuation.

### Alternative: use the maintained container

Replace `run_python` with the [CVMFS container wrapper](index.md#alternative-a-cvmfs-container)
and remove the `conda activate` line. No environment activation is needed
when using the container.
The exact container path is:

```text
/cvmfs/unpacked.cern.ch/ghcr.io/ekatralis/ecloud-containers:latest/
```

Input and output transfers still run in the host worker shell. The environment
recording command will then capture the container's `ECLOUD_CONTAINER_VERSION`.

### Exit codes and checkpoint transfers

The script captures Python's exit code immediately. It uploads the checkpoint
on `0` or `177`, then returns that code. Input or output transfer failures return
`1`. Do not add an unguarded `set -e` around Python: `177` is expected and the
checkpoint must be uploaded before exiting.

The script checks the required files and transfer results, then uploads
`simulation_status.sta` last. Only after all uploads succeed does it delete
the previous `bunch_status_partNN.h5` used to resume the simulation. The new
bunch state and evolution files are kept. If deletion fails, the script logs
a warning and preserves the simulation exit code so continuation can proceed;
the old file can be removed manually.
An interrupted upload can leave files from different parts on EOS. Inspect
the checkpoint before restarting. Run only one workflow per EOS simulation
directory.

## DAGMan continuation

Save as `workflow.dag`:

```text
JOB simulation htcondor.sub
SCRIPT POST simulation post_check.sh $RETURN
RETRY simulation 999 UNLESS-EXIT 0
ABORT-DAG-ON simulation 11 RETURN 1
```

Save as `post_check.sh`:

```bash
#!/bin/bash
case "${1:-missing}" in
    0)   exit 0 ;;
    177) exit 177 ;;
    *)   exit 11 ;;
esac
```

The POST script runs on the submission side. It maps the worker return code
to a DAG decision: `0` completes the node, `177` requests another part, and
any other result becomes `11`, aborting the DAG. The retry count bounds the
number of additional parts; choose it to cover your planned run. See
[DAGMan retry and abort behavior](https://htcondor.readthedocs.io/en/latest/automated-workflows/dagman-advance-functionality.html).

```bash
mkdir -p logs
chmod +x pyparis_executable.sh post_check.sh
condor_submit -dry-run submit.ads htcondor.sub
condor_submit_dag workflow.dag
```

A job held during output transfer can delay the POST script and prevent
continuation. Inspect `condor_q -hold`, the scheduler event log, and the DAGMan
output when continuation stops. This workflow continues **completed parts**;
it does not automatically recover a worker killed mid-part. Resolve incomplete
state or transfer failures before restarting the workflow.

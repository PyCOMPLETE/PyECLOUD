# Buildup simulation scans

Use one independent HTCondor job per scan point. Each job receives common files
and its own input archive, runs in scratch, and returns one explicitly named
archive to EOS. Read the [storage and environment setup](index.md) first.

## Prepare the scan

Keep this submission layout on AFS:

```text
scan/
  htcondor.sub
  job.sh
  sims.txt
  logs/
  point_000/inputs.tgz
  point_001/inputs.tgz
```

`sims.txt` contains one relative directory name per line (`point_000`,
`point_001`, ...), without spaces. Put each point's `.input` and `.beam` files
at the top level of its `inputs.tgz`. Put `run_sim.py`, common geometry, and
other shared files at the top level of `jobfiles.tgz` on EOS. Paths in the
simulation inputs must resolve within the extracted working directory.

## Submit description

Save as `htcondor.sub`. Resource requests below are illustrative; allow scratch
space for inputs, unpacked results, and the compressed archive simultaneously.

```text
universe = vanilla
arguments = ""
executable = job.sh
output = htcondor.out
error = htcondor.err
log = logs/$(ClusterId).log

should_transfer_files = YES
transfer_input_files = $(folder)/inputs.tgz, root://eosproject-e.cern.ch//eos/project/e/ecloud-simulations/YOURDIRECTORY/scan/jobfiles.tgz
when_to_transfer_output = ON_EXIT
transfer_output_files = output.tgz
output_destination = root://eosproject-e.cern.ch//eos/project/e/ecloud-simulations/YOURDIRECTORY/scan/results/$(folder)/
MY.XRDCP_CREATE_DIR = True

+JobFlavour = "tomorrow"
queue folder from sims.txt
```

Each scan point has one EOS directory containing `output.tgz`. Resubmitting
a point reuses that destination. The event log stays on AFS.

## Worker script

Save as `job.sh`. Set `OUTPUT_MAT` to the filename configured by
`filen_main_outp` in your simulation inputs, including the `.mat` extension.
The default option set by PyECLOUD is `Pyecltest.mat`.

```bash
#!/bin/bash
set -e

OUTPUT_MAT="Pyecltest.mat"

# Keep an archive available if a later command fails.
printf 'Job did not complete; see the HTCondor logs.\n' > job-status.txt
tar -czf output.tgz job-status.txt

tar -xzf jobfiles.tgz
rm jobfiles.tgz
tar -xzf inputs.tgz
rm inputs.tgz

# Remove if using container and replace with container command -----
conda activate /eos/user/{yourinitial}/{yourusername}/pyecloud-env
run_python() {
    python "$@"
}
# -------------------------------------------------------------------

run_python -c 'import os, PyECLOUD; key = "ECLOUD_CONTAINER_VERSION"; print(f"{key}={os.environ[key]}" if key in os.environ else f"PyECLOUD={PyECLOUD.__version__}")' > envinfo.txt
run_python run_sim.py

printf 'Simulation completed successfully.\n' > job-status.txt
optional_files=()
if [[ -f stop.txt ]]; then optional_files+=(stop.txt); fi

tar -czvf output.tmp.tgz "$OUTPUT_MAT" progress.txt logfile.txt \
    envinfo.txt job-status.txt "${optional_files[@]}" *.input *.beam
mv output.tmp.tgz output.tgz
```

The simulation runs directly in `_CONDOR_SCRATCH_DIR`. `set -e` stops the
script if extraction, activation, Python, or archiving fails. The final archive
contains the named result and log files, the inputs, and `stop.txt` when present.
Required files must exist for archiving to succeed.

A small fallback archive is created before the simulation. The final archive
replaces it only after `tar` succeeds, so an ordinary failure still leaves
`output.tgz` available for HTCondor to transfer. Forced termination or disk
failure can still prevent delivery. Check the job exit code and logs before
using the results.

### Using the container

Remove `conda activate` and replace `run_python` with the
[CVMFS container wrapper](index.md#alternative-a-cvmfs-container).
The same recording command then writes `ECLOUD_CONTAINER_VERSION` from inside
the container to `envinfo.txt`. With the EOS Conda environment, it records the
installed PyECLOUD version instead. No environment activation is needed for
the container.

## Submit

```bash
mkdir -p logs
condor_submit htcondor.sub
```

The dry run parses the description without submitting jobs; it does not check
worker environment availability or EOS permissions.

## Alternative: transfer individual files from the worker

To avoid archiving, set `transfer_output_files = ""`, remove
`output_destination`, and explicitly upload selected outputs with `xrdcp`
after Python exits. Capture Python's exit code before transferring. Enumerate
existing local files in the shell and check **every** transfer's exit code;
do not pass wildcards to `transfer_output_files`. Fail the job if a required
result is missing or its upload fails. Use the scan point's EOS directory and retain logs for failed simulations too.

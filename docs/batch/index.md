# HTCondor at CERN

Submit buildup scans and single-bunch stability simulations to CERN HTCondor
using the examples below. Replace the example EOS locations,
Conda paths, resources, and simulation inputs before use.

```{toctree}
:maxdepth: 1

buildup-scans
stability
```

## Where files live

For the standard CERN batch setup:

| Location | Use in these examples |
| --- | --- |
| AFS, mounted on submit, scheduler, and worker nodes | Submit descriptions, executable scripts, DAG files, scheduler logs |
| EOS, mounted on submit and worker nodes, but not the standard scheduler nodes | Simulation inputs, results, and optionally the Python environment |
| Worker scratch (`_CONDOR_SCRATCH_DIR`) | Running the simulation and producing local output |

Submit from an AFS directory and keep submission files and scheduler logs there.
Workers can access EOS directly. For HTCondor-managed EOS transfers, use
`root://` URLs with the CERN XRootD transfer plugin, or run `xrdcp` inside
the worker to handle transfers yourself.
See the [CERN batch documentation](https://batchdocs.web.cern.ch/) for more details about HTCondor submissions.

### EOS paths and XRootD endpoints

The examples use `/eos/project/e/ecloud-simulations/YOURDIRECTORY/`.
Replace `YOURDIRECTORY` with your directory within the project.

Choose the XRootD endpoint to match the EOS storage you are using:

| Storage | Endpoint |
| --- | --- |
| EOS user space | `root://eosuser.cern.ch` |
| EOS project space | `root://eosproject-{initial project letter}.cern.ch` |

For the `ecloud-simulations` project, the initial letter is `e`, so the endpoint
is `root://eosproject-e.cern.ch`. A complete file URL looks like:

```text
root://eosproject-e.cern.ch//eos/project/e/ecloud-simulations/YOURDIRECTORY/scan/jobfiles.tgz
```

For user storage, use `root://eosuser.cern.ch` with your `/eos/user/...` path.
Update both the endpoint and the path when switching storage. In the stability
worker script these are `ROOT_URL` and `SIM_PATH`. The script normalizes their
trailing slashes so the combined URL has `//eos/` after the hostname.

## Python environment

Before creating the environment or installing packages, put Conda's package
cache in `/tmp` to avoid filling your AFS home directory:

```bash
mkdir -p /tmp/$USER/conda-pkgs
export CONDA_PKGS_DIRS=/tmp/$USER/conda-pkgs
```

Run these commands in the shell you use for installation.

Create the **complete Conda environment on EOS** so workers can access its
Python interpreter and installed packages. Replace `{yourinitial}` and
`{yourusername}` with your CERN username's first letter and your username:

```bash
conda create -p /eos/user/{yourinitial}/{yourusername}/pyecloud-env python=3.13
conda activate /eos/user/{yourinitial}/{yourusername}/pyecloud-env
```

With this environment active, follow the package installation steps in the
[installation guide](../installation.md), skipping its environment creation
and activation commands.

Batch jobs run in a non-interactive shell where `conda init` has not run, so a bare
`conda activate` fails. The worker scripts therefore load Conda's shell functions first:

```bash
eval "$(conda shell.bash hook)"
conda activate /eos/user/{yourinitial}/{yourusername}/pyecloud-env
```

This requires the `conda` command to be found on the worker. If it is not on the
default `PATH`, call it by its full path, e.g.
`eval "$(/path/to/miniforge3/bin/conda shell.bash hook)"`. Test both lines in a
fresh lxplus session before submitting jobs.

### Alternative: a CVMFS container

We recommend `ghcr.io/ekatralis/ecloud-containers:latest`, maintained for
electron cloud simulations and available through CERN CVMFS. This avoids
maintaining your own Python environment on EOS. The container provides the
Python environment, so no Conda or other environment activation is needed.
Replace environment activation and the Python invocation with a wrapper such as:

```bash
run_python() {
    apptainer exec --cleanenv --env PYTHONNOUSERSITE=1 \
        --home "$_CONDOR_SCRATCH_DIR" --pwd "$PWD" --writable-tmpfs \
        /cvmfs/unpacked.cern.ch/ghcr.io/ekatralis/ecloud-containers:latest/ \
        python "$@"
}
run_python run_sim.py
```

What each option does:

| Option | Effect | Why it is used here |
| --- | --- | --- |
| `--cleanenv` | Starts the container without the host's environment variables. | Host settings such as `PYTHONPATH` or a loaded Conda environment cannot leak in and change which packages Python imports. Remove it only if the job deliberately passes variables into the container. |
| `--env PYTHONNOUSERSITE=1` | Tells Python to ignore packages in `~/.local`. | Packages you installed with `pip install --user` cannot override the container's versions. |
| `--home "$_CONDOR_SCRATCH_DIR"` | Uses the job's scratch directory as the home directory inside the container. | By default Apptainer mounts your AFS home, so caches and configuration files would be written to AFS. With this option they stay in scratch, and the job's files are visible inside the container. |
| `--pwd "$PWD"` | Starts the container in the current directory. | Python runs where the job files were unpacked, so relative input and output paths resolve as they would without the container. |
| `--writable-tmpfs` | Adds a small in-memory writable layer (64 MB by default) on top of the read-only image. | Software that writes caches or temporary files inside the image (for example Matplotlib or Numba caches) would otherwise fail with a read-only file system error. Files written there are lost when the container exits, so results must still go to the working directory. |

The `latest` image evolves; record `ECLOUD_CONTAINER_VERSION` from inside the
container with your run and keep the software environment consistent across
continuation jobs. Run EOS transfers from the host worker shell, where its
credentials and XRootD tools are available.

## Output transfers and held jobs

Both examples save `envinfo.txt`, which records when the simulation started and
the container version or PyECLOUD version it used. The stability example appends
one such record for every simulation part.

Configure output transfers in the **HTCondor submit file** (`htcondor.sub`).
For the `root://` destinations used here, explicitly list individual files in
the submit file's `transfer_output_files`. Do not use wildcards or rely on directory transfer via
this plugin. Archive a directory into one named file, or transfer its files
yourself from the worker.

**Every file listed in `transfer_output_files` must exist when output transfer
runs.** A missing file causes a transfer error and can put the job on hold,
making it appear stuck after the calculation finishes. This is why the buildup
example lists one archive and creates a diagnostic version before starting.
Check the exit code and logs even if the archive was transferred successfully.

With `transfer_output_files = ""`, the worker is responsible for saving the
simulation data. HTCondor still transfers standard output and error to the
locations specified by `output` and `error`.
`when_to_transfer_output = ON_EXIT` does not recover outputs on eviction or
removal. See the [HTCondor file-transfer manual](https://htcondor.readthedocs.io/en/latest/users-manual/file-transfer.html).

Inspect held jobs with `condor_q -hold` and the job's `HoldReason` before
releasing them. Check both the simulation exit status and transfer success.
Start with one short job before submitting a scan or a long continuation chain.

# Running efficient simulations on HTCondor using containers

In this page, we explain an efficient way to run simulations on CERN's computing clusters using containers. We will use the example provided in this repository under:

```text
LHC_6.8TeV_Arcs_Filling_Pattern_Scan/
```

## Table of contents

- [General setup](#general-setup)
- [The Setup](#the-setup)
- [The shared job files](#the-shared-job-files)
- [Generating the simulation inputs](#generating-the-simulation-inputs)
- [The job script](#the-job-script)
- [The submission script](#the-submission-script)
- [Useful commands](#useful-commands)
- [Untarring output files](#untarring-output-files)
- [Transferring files from EOS](#transferring-files-from-eos)

## General setup
Some general technical details which explain why our job submission setup:
- AFS is a (somewhat legacy) shared filesystem, with more limited resources  &rarr; We should avoid using it as much as possible
- EOS is a newer file system designed for large file transfers/accesses &rarr; We should use it to transfer simulation data back and forth
- CVMFS is a software distribution file system designed for frequent accesses &rarr; Using software/containers on CVMFS is ideal
- HTCondor works as a intermediary between the submission nodes (lxplus) and the worker nodes (CERN's clusters that run the jobs). A **Scheduler** controls the negotiation process and matches jobs to workers.
    - A shared FUSE mounted filesystem is a requirement for the jobs to run properly
    - On LxPlus nodes: AFS/EOS/CVMFS are all FUSE mounted
    - On worker nodes: AFS/EOS/CVMFS are all FUSE mounted
    - On Schedulers: Only AFS/CVMFS are FUSE mounted, EOS can only be accessed via the `xrootd` protocol
    - Based on the above, the shared filesystem for all files specified in a submission file **must be** AFS
- HTCondor creates a scratch directory for your jobs to run and can transfer files in and out of this directory as job input/output

> [!NOTE]
> FUSE mounts allow remote filesystems to be mounted and accessed locally as if they were a directory in the local system. 
> For example AFS can be accessed at `/afs/`, EOS at  `/eos/` etc.

## The Setup
Based on the technical details described above our setup must look as follows. We:
- Separate job files that are shared between all jobs (e.g. chamber) and files that change between jobs
- Create a `.tar` archive for all the shared files, which we place on either AFS or EOS depending on the size
- Create our simulation input files on **AFS** and compress them in a `.tar` archive for each job
- Submit from **AFS** and keep logs there
- Compress each job output as a single `.tar` file and transfer that to **EOS**
- Once all jobs finish, run a script that exports all the files from the `.tar` archives and checks for failed simulations

## The shared job files
For the example that we are using, all simulations share the same chamber and photoemission distribution. The script that runs the simulation is also shared between all the simulations. These files can be seen under `LHC_6.8TeV_Arcs_Filling_Pattern_Scan/jobfiles`. 
> [!Important]
> In additions to the files required by the simulation, we also add a file that resembles the `output.tar` that will be transferred on job exit.
> This is important, as missing output files will make the job hang. 
> It is preferrable to let a job finish and determine whether it has failed from the output.
In our case, the `output.mat` file simply contains a `.txt` file that indicates that the job didn't finish properly and to re-run the job. We compress all these files into a single archive called `jobfiles.tar`.

## Generating the simulation inputs
This section will differ depending on the exact submission scripts you have developed, but let's outline the key highlights that are crucial for this step to work regardless of your input generation details. The job generation scripts should:
- Compress the simulation inputs into **a single .tar** file with a **shared name** for all jobs
- Output the paths to all the job inputs as **relative paths** (**not** absolute paths) in a single `.txt` file with each line corresponding to one job.
- Run and create inputs on **AFS only**

> [!Important]
> Ensure that any paths or references to files/inputs are **relative paths** as the directory where the job will run has a randomized name that is not known.
> Specifying an absolute path (e.g. `/afs/...`) means that jobs will run on AFS and **not** in the job scratch directory. (You **WILL** get throttled)

## The job script
The example job file is shown below. Following that, we will explain what each line in the script does:
```bash
#!/bin/bash

CONTAINER_FULLPATH="/cvmfs/unpacked.cern.ch/ghcr.io/ekatralis/ecloud-containers:latest/"
containerrun() {
  apptainer exec \
    --env PYTHONNOUSERSITE=1 \
    --home "$_CONDOR_SCRATCH_DIR" \
    --writable-tmpfs \
    --cleanenv \
    $CONTAINER_FULLPATH \
    "$@"
}

tar -xzf jobfiles.tgz
rm jobfiles.tgz
tar -xzf inputs.tgz
rm inputs.tgz

# Optional: Print node info
echo "************************ NODE INFO *************************"
hostname -A
hostname -I
lscpu
echo "*********************** END NODE INFO ***********************"

# Important: Print container version for future reference
echo "********************** CONTAINER INFO **********************"
containerrun bash -lc 'echo $ECLOUD_CONTAINER_VERSION'
echo "******************** END CONTAINER INFO ********************"

echo "Sim Start --------------------------------------------------"
date

containerrun python run_sim.py 

echo "Sim End ----------------------------------------------------"
date

rm output.tgz
tar -czvf output.tgz output.mat progress.txt logfile.txt $(ls stop.txt 2>/dev/null) *.input *.beam 
```

- We specify the path to the container as an absolute path, as the container is on CVMFS. If we used a custom container that was transferred from EOS, we would use a relative path instead
- The `containerrun` function provides a handy alias for running apptainer with the right arguments:
    - `--env PYTHONNOUSERSITE=1`: Disables python userland, to avoid conflicts with locally installed packages and only use software inside the container
    - `--home "$_CONDOR_SCRATCH_DIR"`: Bind Condor Scratch dir as home inside the container. By default, it would mount your AFS home (**!!**)
    - `--writable-tmpfs`: Provides a 64MB for caches that could end up inside container. By default, container is read-only
    - `--cleanenv`: Prevents environment variable spillover
- We untar the common and job specific archives that were transferred by HTCondor and remove the archives afterwards
- We gather info on the node, optional, but can be useful for troubleshooting problems
- We gather info on the container. **Important**: Preserves reproducibility
- Run the simulation and keep timestamps on start/end times.
- Compress simulation input and output in a single `output.tgz` archive to be transferred back by HTCondor.
    - `$(ls stop.txt 2>/dev/null)` This subshell command is useful for files that could be created by the simulation, but not always

> [!Important]
> Keep all filepaths as **relative paths** unless you explicitly want to access something in the shared filesystem.
> Highly discouraged, unless it is CVMFS.

## The submission script
The example submission file is shown below. Following that, we will once again explain what each line does:
```text
executable              = job.job
arguments               = ""
output                  = htcondor.out
error                   = htcondor.err
log                     = log.$(ClusterId)
should_transfer_files   = YES
transfer_input_files    = $(folder)/inputs.tgz, jobfiles.tgz
output_destination      = root://eosproject-e.cern.ch//eos/project/e/ecloud-simulations/ekatrali/LHC_6.8TeV_Arcs_Filling_Pattern_Scan/$(folder)
when_to_transfer_output = ON_EXIT
transfer_output_files   = output.tgz
MY.XRDCP_CREATE_DIR     = True
+JobFlavour="tomorrow"
queue folder from sims.txt
```
- `executable = job.job`: Select the executable we created in the previous step
- `arguments = ""`: Executable requires no additional arguments
- `output = htcondor.out`: Filename to save `stdout` of the job. This file will end up in EOS together with our output files
- `error = htcondor.err`: Filename to save `stderr` of the job. This file will end up in EOS together with our output files
- `log = log.$(ClusterId)`: Filename to save the log for all jobs. This file can be separate for each job, but it **has to be on AFS**.
- `should_transfer_files = YES`: Ensure that our input files are transferred
- `transfer_input_files = $(folder)/inputs.tgz, jobfiles.tgz`: Transfer inputs for each job as well as shared files:
    - If the shared files are on EOS, the path must be specified using `root://eos{user/project}.cern.ch/.../jobfiles.tgz`
    - In our case the `jobfiles.tgz` is small enough, that the total transfer for all 15000 jobs is around ~2GB, so this file can live on AFS.
- `output_destination = ...`: Set location to transfer the outputs for each job, in this case, we are using an EOS project folder, so the path is set as:
    ```text
    root://eosproject-e.cern.ch//eos/project/e/ecloud-simulations/ekatrali/LHC_6.8TeV_Arcs_Filling_Pattern_Scan/$(folder)
    ```
- `when_to_transfer_output = ON_EXIT`: Transfer output when job finishes
- `MY.XRDCP_CREATE_DIR = True`: Create directories that don't exist on EOS, crucial to avoid errors
- `transfer_output_files = output.tgz`: Specify file that has to be transferred when the job finishes. In this case we transfer the tar archive containing all the simulation data
- `+JobFlavour="tomorrow"`: Set a time limit of 24h for the job
- `queue folder from sims.txt`: Select folders from the relative paths in the `.txt` file that was created by our submission scripts
- Any time `$(folder)` is used, it will be replaced by the values from each line in the `sims.txt` file.

An example of a `sims.txt` file can be seen below:
```bash
$ head sims.txt 
LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoem/LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoemSEY1.00_intens0.30_blength_1.00
LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoem/LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoemSEY1.00_intens0.40_blength_1.00
LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoem/LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoemSEY1.00_intens0.50_blength_1.00
LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoem/LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoemSEY1.00_intens0.60_blength_1.00
LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoem/LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoemSEY1.00_intens0.70_blength_1.00
LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoem/LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoemSEY1.00_intens0.80_blength_1.00
LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoem/LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoemSEY1.00_intens0.90_blength_1.00
LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoem/LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoemSEY1.00_intens1.00_blength_1.00
LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoem/LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoemSEY1.00_intens1.10_blength_1.00
LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoem/LHC_6.8TeV_Arcs_3x48_dipole_Cu2O_pessimistic_photoemSEY1.00_intens1.20_blength_1.00
```

## Useful commands
In these commands replace the log filename with the one that corresponds to the jobs that you ran.

Check the number of unique jobs that finished correctly (should equal to the total number of jobs):
```bash
grep -B1 "return value 0" log.{num} | grep "Job terminated" | awk '{print $2}' | sort -u | wc -l
```
Check for failed jobs:
```bash
grep -E "return value [1-9]|Abnormal termination|Job was evicted|Job was held|Job was aborted" log.{num}
```
Check for jobs that were re-run (usually harmless):
```bash
grep -B1 "return value 0" log.{num} | grep "Job terminated" | awk '{print $2}' | sort | uniq -c | awk '$1>1'
```
Inspect exit codes of specific jobs:
```bash
for j in 10280 10543 11094 5392 7804 13462; do
  echo "=== $j ==="
  grep "11490351.$j" log.{num}
done
```
See full logs of specific jobs:
```bash
for j in 10280 10543 11094 5392; do
  echo "=== $j ==="
  grep -A8 -B2 "005 (11490351.$j.000)" log.{num}
done
```
Get full history and paths of a failed job:
```bash
condor_history 11490351.5607 -long
```
## Untarring output files
Once the simulations have finished running, the archives have to be extracted so that the simulation data can be analyzed. This provides a good opportunity to check for missing output files, unfinished simulations or other problems that could have occured during the job's run. We can do that using the `postjob_untar_and_check.sh` script on a system that has EOS FUSE mounted. For large batches, the script can take a long time to run, as such it is recommended to run it from a `tmux` or `screen` environment to ensure that the connection to that remote node doesn't time out.

Short overview of the script:
- Reads the `sims.txt` file 
- Untars output file in that directory
- Checks that `output.mat` exists
- Checks that the progress file is not below 0.9 (job crashed/finished early)
- If it detects a failed job it prints the filepath in a `failed_sims.txt` file

## Transferring files from EOS
Running large simulation campaigns will often produce very large amounts of data (even hundreds of GBs). While running directly from the EOS fuse mount is a possibility, it is often faster to copy the simulation outputs locally to run the analysis scripts. This can be done via a FUSE mount or SSH connection using rsync/cp or many different methods. The best method to copy files from EOS is actually using the `xrdcp` plugin using parallel streams:
```bash
xrdcp --parallel 8 --retry 3 --verbose -r root://eosproject-e.cern.ch//eos/project/e/ecloud-simulations/{your_path} /local/path/
```
This enables files in directories to be downloaded in parallel, singificantly speeding up the transfer. During my testing, 8 concurrent streams were enough to saturate a gigabit connection. For faster connections, more parallel streams can prove useful, monitor network usage using `btop` and adjust accordingly.
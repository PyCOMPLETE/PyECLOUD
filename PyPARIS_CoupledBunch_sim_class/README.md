# Instructions for running PyECLOUD-PyHEADTAIL coupled bunch simulations on CERN HPC cluster

* Create a workspace folder for software and simulations:
```
mkdir /hpcscratch/user/lmether/sim_workspace
```

* Load MPI: 
```
module load mpi/mvapich2/2.2
```

* Build our software stack in this folder by following the instructions here: https://github.com/PyCOMPLETE/PyECLOUD/wiki/Setup-python-3-%28including-mpi4py%29-without-admin-rights
* Clone this example into the simulation workspace
```
cd /hpcscratch/user/lmether/sim_workspace
git clone https://github.com/lmether/PyECLOUD_coupled_bunch_example
```
* Edit 000s_slurm_job to have the SIM_WORKSPACE variable point to your simulation workspace path

* Launch the simulation:
```
sbatch 000s_slurm_job
```

* Some code to display the simulation results can be found at: https://github.com/giadarol/CoupledBunchAnalysis

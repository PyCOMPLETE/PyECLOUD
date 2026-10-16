#!/usr/bin/bash
#SBATCH -p inf-short
#SBATCH -n 18
#SBATCH -t 1-00:00:00
#SBATCH --hint=nomultithread                                                                                     

module load mpi/mvapich2/2.2

SIM_WORKSPACE=/hpcscratch/user/ecloud/sim_workspace

export PYTHONPATH=$PYTHONPATH:$SIM_WORKSPACE

source $SIM_WORKSPACE/virtualenvs/py2.7/bin/activate

srun python for000_run_sim.py

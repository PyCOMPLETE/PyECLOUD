#!/usr/bin/bash

export PYTHONPATH=$PYTHONPATH:../../../

rm simulation_status.sta

# Run Parallel without MPI
for i in 1 2 3 4
do
    echo "Run $i"
    mpiexec -n 6 --oversubscribe python for000_run_sim.py
done


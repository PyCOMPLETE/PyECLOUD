#!/usr/bin/bash

export PYTHONPATH=$PYTHONPATH:../../../

rm simulation_status.sta

# Run Parallel without MPI
for i in 1 2 3 4
do
    echo "Run $i"
    ../../../PyPARIS/multiprocexec.py -n 6 sim_class=PyPARIS_CoupledBunch_sim_class.Simulation.Simulation --multiturn
done

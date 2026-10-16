#!/usr/bin/bash

export PYTHONPATH=$PYTHONPATH:../../../

rm simulation_status.sta

# Run Parallel without MPI
for i in 1 2 3
do
    echo "Run $i"
    python -m PyPARIS.multiprocexec -n 3 "sim_class=PyPARIS_sim_class.Simulation.Simulation(param_file='param_folder/Simulation_parameters.py')"
done


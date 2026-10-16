#!/usr/bin/bash

export PYTHONPATH=$PYTHONPATH:../../../

rm simulation_status.sta

# Run MPI
for i in 1 2 3
do
    echo "Run $i"
    mpiexec -n 3 python -m PyPARIS.withmpi sim_class=PyPARIS_sim_class.Simulation.Simulation
done

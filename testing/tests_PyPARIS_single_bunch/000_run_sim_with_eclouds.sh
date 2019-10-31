#!/usr/bin/bash

export PYTHONPATH=../../../PyPARIS/:$PYTHONPATH

# Run Parallel without MPI
# ../../../PyPARIS/multiprocexec.py -n 3 sim_class=Simulation_with_eclouds.Simulation

# Run Serial
#../../../PyPARIS/serialexec.py sim_class=Simulation_with_eclouds.Simulation

# Run MPI
mpiexec -n 4 ../../../PyPARIS/withmpi.py sim_class=Simulation_with_eclouds.Simulation

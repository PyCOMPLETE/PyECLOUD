#!/usr/bin/bash

# Run MPI
# mpiexec -n 16 python for000_run_sim.py

# Run multiproc
../../../PyPARIS/multiprocexec.py -n 16 sim_class=Simulation.Simulation --multiturn

# Run single core
#python for000_run_sim.py

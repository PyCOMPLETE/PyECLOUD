#!/usr/bin/bash

export PYTHONPATH=$PYTHONPATH:../../../

# Run Serial
python -m PyPARIS.serialexec sim_class=PyPARIS_sim_class.Simulation.Simulation

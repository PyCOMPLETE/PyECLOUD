#!/usr/bin/bash

set -euo pipefail

# Run Parallel without MPI
python -m PyPARIS.multiprocexec -n 3 sim_class=Simulation_with_eclouds.Simulation

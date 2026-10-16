#!/usr/bin/bash

set -euo pipefail

# Run Serial
python -m PyPARIS.serialexec sim_class=Simulation_with_eclouds.Simulation

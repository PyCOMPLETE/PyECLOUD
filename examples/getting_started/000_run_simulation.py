"""Run a simulation from a folder containing PyECLOUD input files."""

import argparse
from PyECLOUD.buildup_simulation import BuildupSimulation

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('input_folder', nargs='?', default='.')
args = parser.parse_args()

sim = BuildupSimulation(pyecl_input_folder=args.input_folder)
sim.run()

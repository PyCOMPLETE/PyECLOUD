"""Resume a saved simulation state using its original input folder."""

import argparse
from PyECLOUD.buildup_simulation import BuildupSimulation

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('input_folder')
parser.add_argument('state_file', help='Path to the saved simulation_state_*.pkl file')
args = parser.parse_args()

sim = BuildupSimulation(pyecl_input_folder=args.input_folder)
sim.load_state(args.state_file, load_from_folder='')
sim.run()

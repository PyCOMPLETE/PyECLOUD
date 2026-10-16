from PyECLOUD.buildup_simulation import BuildupSimulation
import argparse
from os import path

parser = argparse.ArgumentParser()
parser.add_argument('--input_folder',type=str, default="./" , help='Folder containing input files')
parser.add_argument('--output_folder',type=str, help='Folder containing output file')
args = parser.parse_args()

if args.output_folder:
    out_file = path.join(args.output_folder, "output.mat")
else:
    out_file = path.join(args.input_folder, "output.mat")

sim = BuildupSimulation(pyecl_input_folder=args.input_folder, filen_main_outp=out_file)
sim.run()

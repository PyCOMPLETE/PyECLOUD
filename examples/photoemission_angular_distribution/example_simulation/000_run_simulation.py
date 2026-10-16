
import sys
import os
import time
import argparse
BIN = os.path.expanduser("../../../../")  # folder containing PyECLOUD, PyPIC, PyKLU
if BIN not in sys.path:
    sys.path.append(BIN)

from PyECLOUD.buildup_simulation import BuildupSimulation


time_0 = time.time()
sim = BuildupSimulation()
sim.run()

time_needed = time.time() - time_0


print('')
print('Simulation done in %.2f s!' % time_needed)
print('')

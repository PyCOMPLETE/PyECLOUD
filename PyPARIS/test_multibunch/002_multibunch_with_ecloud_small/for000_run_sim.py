# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

import os, sys
import sys, os
BIN = os.path.expanduser("../../../")
sys.path.append(BIN)


from PyPARIS.ring_of_CPUs_multiturn import RingOfCPUs_multiturn
import Simulation 

simulation_content = Simulation.Simulation()
myCPUring = RingOfCPUs_multiturn(simulation_content)
myCPUring.run()

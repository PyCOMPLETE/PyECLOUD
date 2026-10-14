# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

from . import Simulation as sim_mod

ring = sim_mod.get_serial_CPUring()
ring.run()

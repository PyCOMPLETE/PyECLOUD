# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

#!/usr/bin/env python
import multiprocessing as mp
import numpy as np
import os, sys
import importlib

from . import parse_sim_class_string as psc

sys.path.append('./')

sim_module_string = sys.argv[1].split('=', 1)[1]

module_name, class_name, dict_kwargs = psc.parse_sim_class_string(
        sim_module_string)

SimModule = importlib.import_module(module_name)
SimClass = getattr(SimModule, class_name)
simulation_content = SimClass(**dict_kwargs)

from PyPARIS.ring_of_CPUs import RingOfCPUs
myCPUring = RingOfCPUs(simulation_content)
myCPUring.run()

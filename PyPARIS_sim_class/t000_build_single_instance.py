# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

from . import Simulation as sim_mod

N_proc_pretend = 4
myid_pretend = 2

sim = sim_mod.get_sim_instance(N_proc_pretend, myid_pretend)


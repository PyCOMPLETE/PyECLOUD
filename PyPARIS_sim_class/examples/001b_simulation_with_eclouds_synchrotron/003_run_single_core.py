import sys
sys.path.append('../../../')

import PyPARIS_sim_class.Simulation as sim_mod

ring = sim_mod.get_serial_CPUring()
ring.run()

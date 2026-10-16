import sys
sys.path.append('../../../')

import PyPARIS_sim_class.Simulation as sim_mod
import PyPARIS.util as pu

sim_content = sim_mod.Simulation(
    param_file='param_folder/Simulation_parameters.py')

ring = pu.get_serial_CPUring(sim_content, init_sim_objects_auto=True)
ring.run()

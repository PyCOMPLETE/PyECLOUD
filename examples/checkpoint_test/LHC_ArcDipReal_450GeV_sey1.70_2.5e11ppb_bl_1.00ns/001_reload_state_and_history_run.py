import sys
import os
BIN = os.path.expanduser("../../../../")
sys.path.append(BIN)

from PyECLOUD.buildup_simulation import BuildupSimulation

saved_states_list = [f for f in os.listdir('./') if (f.find('simulation_checkpoint_') != -1)]
if len(saved_states_list) == 0:
    raise Exception('No simulation checkpoint file was found')
saved_states_list.sort()


last_saved_state = saved_states_list[-1]

sim = BuildupSimulation()
sim.load_checkpoint(last_saved_state)
sim.run()
# merge_Pyecltest.merge_Pyecltest()

import sys, os
BIN = os.path.expanduser("../")
sys.path.append(BIN)


from PyECLOUD.buildup_simulation import BuildupSimulation


sim = BuildupSimulation()
sim.load_state('simulation_state_0.pkl')
sim.run()

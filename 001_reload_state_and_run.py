import sys, os
BIN = os.path.expanduser("../")
sys.path.append(BIN)


from PyECLOUD.buildup_simulation import BuildupSimulation


sim = BuildupSimulation()

# Optionally enable video saving
sim.pyeclsaver.flag_video = False
sim.pyeclsaver.flag_sc_video = False

sim.load_state('simulation_state_0.pkl')
sim.run()

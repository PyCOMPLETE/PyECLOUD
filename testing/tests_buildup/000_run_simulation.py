import sys, os
BIN = os.path.expanduser("../../../")
sys.path.append(BIN)

sim_folder = 'LHC_ArcDipReal_450GeV_sey1.60_2.5e11ppb_bl_1.00ns'
# sim_folder = 'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns'


from PyECLOUD.buildup_simulation import BuildupSimulation


sim = BuildupSimulation(pyecl_input_folder = sim_folder)
sim.run()

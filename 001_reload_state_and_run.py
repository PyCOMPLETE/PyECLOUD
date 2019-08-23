import sys, os
BIN = os.path.expanduser("../")
sys.path.append(BIN)


from PyECLOUD.buildup_simulation import BuildupSimulation
input_folder = 'testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_stress_saver/'

sim = BuildupSimulation(pyecl_input_folder=input_folder)

# Optionally enable video saving
for cloud in sim.cloud_list:
    cloud.pyeclsaver.flag_video = False
    cloud.pyeclsaver.flag_sc_video = False

sim.load_state(input_folder+'/simulation_state_0.pkl')
sim.run()

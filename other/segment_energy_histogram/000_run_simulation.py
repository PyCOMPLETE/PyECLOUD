import numpy as np

from PyECLOUD.buildup_simulation import BuildupSimulation

sim_input_folder = '../../testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns'

sim = BuildupSimulation(pyecl_input_folder=sim_input_folder, 
        filen_main_outp='./Pyecltest.mat',
        extract_sey=False)
sim.run(t_end_sim=25e-9)

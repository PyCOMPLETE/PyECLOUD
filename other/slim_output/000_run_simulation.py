import numpy as np

from PyECLOUD.buildup_simulation import BuildupSimulation

sim_input_folder = '../../testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns'

save_only = 't Nel_timep xg_hist nel_hist'.split()

sim = BuildupSimulation(pyecl_input_folder=sim_input_folder, 
        filen_main_outp='./Pyecltest.mat',
        flag_En_hist_seg = True, 
        extract_sey=False,
        Nbin_En_hist= 300, 
        En_hist_max= 2500.,  #eV i
        save_only = save_only
        )

sim.run(t_end_sim=None)

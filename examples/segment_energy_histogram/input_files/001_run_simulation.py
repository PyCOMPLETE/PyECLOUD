import numpy as np

from PyECLOUD.buildup_simulation import BuildupSimulation

sim_input_folder = '../../testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns'

save_once_custom_observables = {
        'Vx': lambda sim: sim.cloud_list[0].impact_man.chamb.Vx,
        'Vy': lambda sim: sim.cloud_list[0].impact_man.chamb.Vy,
        'L_edg': lambda sim: sim.cloud_list[0].impact_man.chamb.L_edg,
        }


sim = BuildupSimulation(pyecl_input_folder=sim_input_folder, 
        filen_main_outp='./Pyecltest.mat',
        filename_chm='chamber.mat',
        flag_En_hist_seg = True, 
        extract_sey=False,
        Nbin_En_hist= 300, 
        En_hist_max= 6000.,  #eV 
        save_once_custom_observables=save_once_custom_observables)

sim.run(t_end_sim=300e-9)

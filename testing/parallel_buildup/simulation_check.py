import numpy as np

from PyECLOUD.buildup_simulation import BuildupSimulation

sim_input_folder = '../../testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns'


from mpi4py import MPI
comm = MPI.COMM_WORLD

myid = comm.Get_rank()

x_aper = 2.3e-2
x_max_list = [-1e-2, 0, 1e-2, x_aper]

sim = BuildupSimulation(pyecl_input_folder=sim_input_folder,
        filen_main_outp='./Pyecltest_%02d.mat'%myid,
        flag_En_hist_seg = True,
        extract_sey=False,
        Nbin_En_hist= 300,
        En_hist_max= 2500.,  #eoV
        x_max_init_unif = x_max_list[myid]
        )
sim.spacech_ele.comm = comm

sim.run(t_end_sim=5*25e-9)

ec = sim.cloud_list[0]

import matplotlib.pyplot as plt
plt.close('all')
fig = plt.figure(100 + myid)
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)

ax1.plot(ec.pyeclsaver.xg_hist, np.sum(ec.pyeclsaver.nel_hist, axis=0))
ax2.plot(sim.spacech_ele.xg, np.sum(sim.spacech_ele.rho, axis=1))


plt.show()

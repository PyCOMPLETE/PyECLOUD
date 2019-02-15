from PyECLOUD.buildup_simulation import BuildupSimulation


sim = BuildupSimulation()

t_end = 5e-9
sim.run(t_end_sim = t_end)

'''
ec = sim.cloud_list[0]
plot(ec.MP_e.x_mp[:ec.MP_e.N_mp], ec.MP_e.y_mp[:ec.MP_e.N_mp], '.')
'''

import numpy as np
import matplotlib.pyplot as plt
from PyECLOUD.dynamics_Boris_f2py import pusher_Boris
from scipy.constants import e, m_e

pusher = pusher_Boris(Dt=1e-6, B0x=0, B0y=0, B0z=0, B_map_file="Bfile_onlyBz.mat", fact_Bmap=1)

MP_e = type('test', (), {})()
MP_e.N_mp = 3
MP_e.x_mp = np.linspace(0, 0.1, MP_e.N_mp)
MP_e.y_mp = np.zeros_like(MP_e.x_mp)
MP_e.z_mp = np.zeros_like(MP_e.x_mp)
MP_e.vx_mp = np.zeros_like(MP_e.x_mp)
MP_e.vy_mp = np.zeros_like(MP_e.x_mp) + 10000
MP_e.vz_mp = np.zeros_like(MP_e.x_mp)
MP_e.charge = -e
MP_e.mass = m_e

Ex_n = np.zeros_like(MP_e.x_mp)
Ey_n = np.zeros_like(MP_e.x_mp)

plt.figure(1)
for ii in range(100):
    for jj in range(100):
        pusher.step(MP_e, Ex_n, Ey_n)
    plt.plot(MP_e.x_mp, MP_e.y_mp, 'k.')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
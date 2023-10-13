import numpy as np
import matplotlib.pyplot as plt
from PyECLOUD.dynamics_Boris_f2py import B_file

myfile = B_file(B0x=0., B0y=0., B0z=0., fact_Bmap=1., B_map_file="Bfile.mat")

xobs = 0.1
yobs = 0.1

xp = np.linspace(-6, 6, 100)
Bx_x = np.zeros_like(xp)
By_x = np.zeros_like(xp)
Bz_x = np.zeros_like(xp)
for ii in range(len(xp)):
    Bx_x[ii], By_x[ii], Bz_x[ii] = myfile.get_B(xp[ii], yobs)


yp = np.linspace(-6, 6, 100)
Bx_y = np.zeros_like(yp)
By_y = np.zeros_like(yp)
Bz_y = np.zeros_like(yp)
for ii in range(len(yp)):
    Bx_y[ii], By_y[ii], Bz_y[ii] = myfile.get_B(xobs, yp[ii])

plt.figure(1)
plt.plot(xp, Bx_x,'.', label='Bx')
plt.plot(xp, By_x,'.', label='By')
plt.plot(xp, 1e4*Bz_x,'.', label='$10^{4}$ Bz')
plt.xlabel('x')
plt.legend()

plt.figure(2)
plt.plot(yp, Bx_y,'.', label='Bx')
plt.plot(yp, By_y,'.', label='By')
plt.plot(yp, 1e4*Bz_y,'.', label='$10^{4}$ Bz')
plt.xlabel('y')
plt.legend()
plt.show()
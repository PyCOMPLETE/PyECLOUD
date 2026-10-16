import os
import numpy as np
import matplotlib.pyplot as plt

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), 'comparison_plots')
os.makedirs(OUTPUT_DIR, exist_ok=True)
from PyECLOUD.dynamics_Boris_f2py import B_file

myfile = B_file(B0x=0., B0y=0., B0z=0., fact_Bmap=1., B_map_file="Bfile.mat")

xobs = 0.1
yobs = 0.1

xp = np.linspace(-6, 6, 100)
Bx_x = np.zeros_like(xp)
By_x = np.zeros_like(xp)
Bz_x = np.zeros_like(xp)
for ii in range(len(xp)):
    Bx, By, Bz = myfile.get_B(xp[ii], yobs)
    Bx_x[ii], By_x[ii], Bz_x[ii] = (
        np.asarray(v).item() for v in (Bx, By, Bz)
    )


yp = np.linspace(-6, 6, 100)
Bx_y = np.zeros_like(yp)
By_y = np.zeros_like(yp)
Bz_y = np.zeros_like(yp)
for ii in range(len(yp)):
    Bx, By, Bz = myfile.get_B(xobs, yp[ii])
    Bx_y[ii], By_y[ii], Bz_y[ii] = (
        np.asarray(v).item() for v in (Bx, By, Bz)
    )

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
plt.figure(1).savefig(os.path.join(OUTPUT_DIR, 'figure_01.png'), dpi=200)
plt.figure(2).savefig(os.path.join(OUTPUT_DIR, 'figure_02.png'), dpi=200)
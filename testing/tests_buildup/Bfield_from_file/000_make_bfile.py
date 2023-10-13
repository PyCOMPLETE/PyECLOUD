import numpy as np
import scipy.io as sio


def Bx_f(x,y):
    return 6 * x**2 * y**2 

def By_f(x,y):
    return 0.1 * x * y 

def Bz_f(x,y):
    return 1e-6 * (1 + 0.1*x**4)


xmax = 5.0
ymax = 5.0
Nx = 100
Ny = 100
xx = np.linspace(-xmax, xmax, Nx)
yy = np.linspace(-ymax, ymax, Ny)
Bx = np.zeros([Nx, Ny])
By = np.zeros([Nx, Ny])
Bz = np.zeros([Nx, Ny])

for ix in range(Nx):
    for iy in range(Ny):
        Bx[ix, iy] = Bx_f(xx[ix], yy[iy])
        By[ix, iy] = By_f(xx[ix], yy[iy])
        Bz[ix, iy] = Bz_f(xx[ix], yy[iy])

B_dict = {"xx" : xx, "yy" : yy,
          "Bx" : Bx, "By" : By,
          "Bz" : Bz}
B_dict_onlyBz = {"xx" : xx, "yy" : yy,
          "Bx" : 0*Bx, "By" : 0*By,
          "Bz" : Bz}
B_dict_noBz = {"xx" : xx, "yy" : yy,
                "Bx" : Bx, "By" : By,
                }

sio.savemat("Bfile.mat", B_dict)
sio.savemat("Bfile_noBz.mat", B_dict_noBz)
sio.savemat("Bfile_onlyBz.mat", B_dict_onlyBz)
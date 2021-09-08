import numpy as np
import matplotlib.pyplot as plt

E_impact_eV = np.linspace(0, 1000., 10000)
theta = 0 * np.pi/180 #0.

Ehat_e = 0
W = 100
p = 0.9
P1e_inf = 0.07
P1e_hat = 0.5
e1 = 0.26
e2 = 2

P1r_inf = 0.74
E_r = 40.
r = 1
r1 = 0.26
r2 = 2

delta_e_0 = P1e_inf + (P1e_hat-P1e_inf)*np.exp(-((np.abs(E_impact_eV-Ehat_e)/W)**p)/p)
delta_r_0 = P1r_inf * (1. - np.exp(-(E_impact_eV/E_r)**r))

delta_e = delta_e_0 * (1. + e1*(1-(np.cos(theta))**e2))
delta_r = delta_r_0 * (1. + r1*(1-(np.cos(theta))**r2))

plt.close('all')
fig1 = plt.figure(1)
ax = fig1.add_subplot(111)
ax.plot(E_impact_eV, delta_e)
ax.plot(E_impact_eV, delta_r)
ax.plot(E_impact_eV, delta_e + delta_r)

ax.grid(True)

plt.show()

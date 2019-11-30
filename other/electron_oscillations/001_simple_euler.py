import numpy as np

N_steps = 100
k_field = -0.05

Dt = 1.

x =1e-3
vx = 0.

x_record = []

for ii in range(N_steps):
    vx += k_field*x*Dt
    x += vx*Dt
    x_record.append(x)

import matplotlib.pyplot as plt
plt.close('all')
plt.plot(x_record)
plt.grid(True)
plt.show()


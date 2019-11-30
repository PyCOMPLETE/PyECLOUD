import numpy as np

N_steps = 100
N_substeps = 10
k_field = -0.05

Dt = 1.

x =1e-3
vx = 0.

x_record = []

for ii in range(N_steps):
    E = k_field*x
    for ssn in range(N_substeps):
        vx += E*Dt/N_substeps 
        x += vx*Dt/N_substeps
    x_record.append(x)

import matplotlib.pyplot as plt
plt.close('all')
plt.plot(x_record)
plt.grid(True)
plt.show()


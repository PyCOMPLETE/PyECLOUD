import numpy as np
import matplotlib.pyplot as plt

import PyECLOUD.myfilemanager as mfm

i_det = 6

ob = mfm.myloadmat_to_obj('Pyecltest.mat')

plt.close('all')
fig1 = plt.figure(1)
sp1 = fig1.add_subplot(1,1,1)
sp1.plot(ob.t/ob.b_spac, ob.Nel_timep)

fig2 = plt.figure()
sp21 = fig2.add_subplot(2,1,1)
mask_det = np.logical_and(ob.t>=i_det*ob.b_spac, ob.t<(i_det+1)*ob.b_spac) 
t_det = ob.t[mask_det]
t_det -= t_det[0]

sp21.plot(t_det, ob.Nel_timep[mask_det]/ob.Nel_timep[mask_det][0])


# Compute simple model
N_b = 1.15e11
R = 2.e-2

from scipy.constants import e as qe
from scipy.constants import c as clight
from scipy.constants import epsilon_0
from scipy.constants import m_e
Dp_bar = 3./4. * qe**2 * N_b/(np.pi*epsilon_0*clight*R)
E_bar_J = Dp_bar**2/(2*m_e)
E_bar_eV = E_bar_J/qe

import PyECLOUD.sec_emission_model_ECLOUD as seymod

del_bar, _ = seymod.yield_fun2(
        E=np.atleast_1d(E_bar_eV), costheta=1., Emax=ob.Emax, 
        del_max=ob.del_max, R0=ob.R0, E0=ob.E0, s=ob.s, 
        flag_costheta_delta_scale=True, flag_costheta_Emax_shift=True)

def phi_emit(E):
    return 1./E*1./(ob.sigmafit*np.sqrt(2*np.pi))*np.exp(-(np.log(E)-ob.mufit)**2/(2*ob.sigmafit**2))

k_ele = qe**2/(2*np.pi*epsilon_0*R**2*m_e)

N_vect = np.logspace(8, 11, 100)

N = N_vect[-1]

tanh2 = (np.tanh(0.5*np.sqrt(k_ele*N)*ob.b_spac))**2

E_minus_eV = 0.5*m_e*k_ele*N*R**2*tanh2/qe
E_plus_eV = 0.5*m_e*k_ele*N*R**2/tanh2/qe

E_vect = np.linspace(0, 50, 1000)

fig100 = plt.figure(100)
sp100 = fig100.add_subplot(1,1,1)
sp100.plot(E_vect, phi_emit(E_vect))

plt.show()

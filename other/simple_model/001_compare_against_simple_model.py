import numpy as np
import matplotlib.pyplot as plt

import PyECLOUD.myfilemanager as mfm
import PyECLOUD.mystyle as ms
i_det = 6

ob = mfm.myloadmat_to_obj('Pyecltest.mat')

plt.close('all')
ms.mystyle_arial(fontsz=14, dist_tick_lab=5, traditional_look=False)
fig1 = plt.figure(1)
sp1 = fig1.add_subplot(1,1,1)
sp1.semilogy(ob.t/ob.b_spac, ob.Nel_timep)

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

# Compute representative kick
Dp_bar = 3./4. * qe**2 * N_b/(np.pi*epsilon_0*clight*R)*2/3
E_bar_J = Dp_bar**2/(2*m_e)
E_bar_eV = E_bar_J/qe

import PyECLOUD.sec_emission_model_ECLOUD as seymod

# Compute corresponding SEY
del_bar, _ = seymod.yield_fun2(
        E=np.atleast_1d(E_bar_eV), costheta=1., Emax=ob.Emax, 
        del_max=ob.del_max, R0=ob.R0, E0=ob.E0, s=ob.s, 
        flag_costheta_delta_scale=True, flag_costheta_Emax_shift=True)
del_bar = del_bar[0]
# Compute electric field gradient coefficient
k_ele = qe**2/(2*np.pi*epsilon_0*R**2*m_e)

# Compute energy spectrum
def phi_emit(E):
    return 1./E*1./(ob.sigmafit*np.sqrt(2*np.pi))*np.exp(-(np.log(E)-ob.mufit)**2/(2*ob.sigmafit**2))


def compute_P_surv(N):
    tanh2 = (np.tanh(0.5*np.sqrt(k_ele*N)*ob.b_spac))**2

    E_minus_eV = 0.5*m_e*k_ele*N*R**2*tanh2/qe
    E_plus_eV = 0.5*m_e*k_ele*N*R**2/tanh2/qe

    E_vect = np.linspace(E_minus_eV, E_plus_eV, 1000)
    P_surv = np.trapz(phi_emit(E_vect), E_vect)

    return P_surv, E_minus_eV, E_plus_eV


# Evaluate abosorption factor for different electron densities
N_vect = np.logspace(6, 11, 200)
P_surv_vect = N_vect*0
E_minus_vect = N_vect*0
E_plus_vect = N_vect*0
for ii, N in enumerate(N_vect):

    P_surv, E_minus_eV, E_plus_eV = compute_P_surv(N)

    P_surv_vect[ii] = P_surv
    E_minus_vect[ii] = E_minus_eV
    E_plus_vect[ii] = E_plus_eV


# Identify saturation level
N_sat = N_vect[np.argmin(np.abs(1./del_bar-P_surv_vect))]

# Buildup from the model
N_pass = np.int(np.max(ob.t/ob.b_spac))

N_0 = ob.Nel_timep[0]
N_iter_vect = [N_0]
del_eff_vect = []
N_after_impact = []
P_surv_iter_vect = []
for ii in range(N_pass):
    P_surv, _, _ = compute_P_surv(N_iter_vect[-1])
    P_surv_iter_vect.append(P_surv)
    del_eff_vect.append(del_bar*P_surv)
    N_after_impact.append(N_iter_vect[-1]*del_bar)
    N_iter_vect.append(N_iter_vect[-1]*del_bar*P_surv)
# Some plots
fig3 = plt.figure(3)
fig3.set_facecolor('w')
sp31 = fig3.add_subplot(2,1,1)
sp32 = fig3.add_subplot(2,1,2, sharex=sp31)
sp31.semilogx(N_vect, 1/P_surv_vect)
sp31.axhline(del_bar, linestyle='--', linewidth=2, color='k')

sp32.semilogx(N_vect, E_minus_vect)
sp32.semilogx(N_vect, E_plus_vect)

sp31.set_ylim(1, 2.)
sp32.set_ylim(0, 10)

for sp in [sp31, sp32]:
    sp.grid(True)

sp31.axvline(N_sat)

sp1.axhline(N_sat)
sp1.plot(N_iter_vect, '.')
sp1.plot(N_after_impact, '.r')
# ## A check
# E_vect = np.linspace(0, 50, 1000)
# 
# fig100 = plt.figure(100)
# sp100 = fig100.add_subplot(1,1,1)
# sp100.plot(E_vect, phi_emit(E_vect))

# Figure for Handbook
fig4 = plt.figure(4)
sp1 = plt.subplot2grid(fig=fig4, shape=(1,3), loc=(0,0), colspan=2)
sp2 = plt.subplot2grid(fig=fig4, shape=(1,3), loc=(0,2), colspan=1,
                       sharey=sp1)

sp1.semilogx(N_vect, E_minus_vect, color='C0')
sp1.semilogx(N_vect, E_plus_vect, color='C0')
sp1.fill_between(x=N_vect, y1=E_minus_vect, y2=E_plus_vect,
                 alpha=0.5)

E = np.linspace(0, 25, 1000)
sp2.plot(phi_emit(E), E)
sp2.fill_betweenx(x2=phi_emit(E), x1=0, y=E, alpha=0.5)

sp1.set_xlim(1e6, 3e10)
sp1.set_ylim(0, 20)
sp2.set_xlim(left=0)
sp1.set_ylabel('Electron energy (eV)')
sp1.set_xlabel(r'Electron  density (m$^{-1}$)')
sp2.set_xlabel(r'$\phi_{emit}(E)$')
plt.setp(sp2.get_yticklabels(), visible=False)
fig4.subplots_adjust(bottom=.265, wspace=.063)

fig5 = plt.figure()
ax5 = fig5.add_subplot(111)
i_vect = np.arange(0, 100)
ax5.semilogy(ob.t/ob.b_spac, ob.Nel_timep, color='grey', alpha=.7,
             label='simulation')
ax5.semilogy(i_vect, 3*N_vect[0]*del_eff_vect[0]**i_vect,
             linestyle='--', linewidth=2, color='C3',
             label=r'$N_0\delta_{eff}^i$')
ax5.axhline(y=N_sat,
            linestyle='--', linewidth=2, color='C0', label=r'$N_{sat}$')
ax5.set_xlim(0, 80)
ax5.set_ylim(1e6, 1e10)
ax5.set_ylabel(r'Electron  density (m$^{-1}$)')
ax5.set_xlabel('Bunch passage')
ax5.legend(loc='lower right', frameon=False)
fig5.subplots_adjust(bottom=.265, left=.15)



plt.show()

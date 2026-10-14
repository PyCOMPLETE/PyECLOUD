import numpy as np
import matplotlib.pyplot as plt

from PyECLOUD import geom_impact_poly_fast_impact as gipfi
import PyECLOUD.sec_emission_model_ECLOUD_nunif as se
import PyECLOUD.mystyle as ms

chamb = gipfi.polyg_cham_geom_object('chamber.mat', flag_non_unif_sey=True)

semod = se.SEY_model_ECLOUD_non_unif_charging(chamb, 
        Emax=-1., del_max=-1., R0=-1, E0=150.,
        E_th=35., sigmafit=1.0828, mufit=1.6636,
        switch_no_increase_energy=0, 
        thresh_low_energy=-1, 
        secondary_angle_distribution='cosine_3D')

Q_Qmax_test_vect = np.linspace(0, 1, 6)

E_test = np.linspace(0, 1600, 10000)

plt.close('all')
ms.mystyle_arial(fontsz=16)

fig1 = plt.figure(1)
fig1.set_facecolor('w')
ax1 = fig1.add_subplot(1,1,1)

for iq, Q_Qmax in enumerate(Q_Qmax_test_vect):

    colorcurr = ms.colorprog(iq, len(Q_Qmax_test_vect), cm='rainbow')
   
    semod.Q_segments[:] = Q_Qmax*semod.Q_max_segments
    
    yiel, flag_elast, flag_truesec = semod.SEY_process(
            nel_impact=E_test*0+1.,
            E_impact_eV=E_test, 
            costheta_impact=E_test*0+1.,
            i_impact=np.int_(E_test*0+1))
    
    ax1.plot(E_test, yiel, label='%.2f'%Q_Qmax, linewidth=2,
                color=colorcurr)

ax1.grid(True)
ax1.set_xlabel('Energy [eV]')
ax1.set_ylabel('Secondary Electron Yield')
ax1.legend(loc='lower right', prop={'size':16}, title='Q/Qmax',
        ncol=2)

fig1.subplots_adjust(
        bottom=.14
        )
# Check charge accumulation
semod.Q_segments[:] = 0.

plt.figure(2)

for ii in range(20):

    yiel, flag_elast, flag_truesec = semod.SEY_process(nel_impact=E_test*0+1.,
            E_impact_eV=E_test, 
            costheta_impact=E_test*0+1.,
            i_impact=np.int_(E_test*0+1))
    
    _, _, _ = semod.SEY_process(nel_impact=5e9*np.ones(1000)/1000.,
            E_impact_eV=300*np.ones(1000),
            costheta_impact=1.*np.ones(1000),
            i_impact=np.int_(np.ones(1000)))
    
    plt.plot(E_test, yiel)


plt.show()


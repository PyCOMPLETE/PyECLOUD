# generate example figures with main PyECLOUD outputs
# Thanks for E. Belli, L. Mether and A. Romano for preparing this script

import sys, os
BIN = os.path.expanduser("../../../")
sys.path.append(BIN)

import pylab as pl
import numpy as np

import PyECLOUD.myloadmat_to_obj as mlo
import PyECLOUD.mystyle as ms

pl.close('all')
ms.mystyle_arial(fontsz=16)
dpiset = 200

ob=mlo.myloadmat_to_obj('FCC_Quad_25ns_50.00TeV_1e-04_R0.1_sey1.5/Pyecltest_ref.mat')

ifig = 0

#8.
ifig+=1; pl.figure(ifig)
pl.plot(ob.t, ob.lam_t_array)
pl.xlabel('Time [s]')
pl.ylabel('Beam profile [p/m]') 
ms.scix(); pl.grid('on')
pl.suptitle('Var. name: lam_t_array\nBeam density at each time step')
pl.subplots_adjust(top=.82, bottom=.14)
pl.savefig('fig%02d.png'%ifig, dpi=dpiset)

#6.
ifig+=1; pl.figure(ifig)
pl.plot(ob.t, ob.Nel_timep)
pl.xlabel('Time [s]')
pl.ylabel('Number of $e^-$ per unit length [$m^{-1}$]') 
ms.scix(); pl.grid('on')
pl.suptitle('Var. name: Nel_timep\nNumber of electrons in the chamber at each time step')
pl.subplots_adjust(top=.82, bottom=.14)
pl.savefig('fig%02d.png'%ifig, dpi=dpiset)

#5.
ifig+=1; pl.figure(ifig)
pl.plot(ob.t, ob.Nel_imp_time)
pl.xlabel('Time [s]')
pl.ylabel('Number of impacting $e^-$ per unit length [$m^{-1}$]') 
ms.scix(); pl.grid('on')
pl.suptitle('Var. name: Nel_imp_time\nNumber of electrons that impact the walls in each time-step')
pl.subplots_adjust(top=.82, bottom=.14)
pl.savefig('fig%02d.png'%ifig, dpi=dpiset)

#4.
ifig+=1; pl.figure(ifig)
pl.plot(ob.t, ob.Nel_emit_time)
pl.xlabel('Time [s]')
pl.ylabel('Number of emitted $e^-$ per unit length [$m^{-1}$]') 
ms.scix(); pl.grid('on')
pl.suptitle('Var. name: Nel_emit_time\nNumber of electrons emitted by walls in each time-step')
pl.subplots_adjust(top=.82, bottom=.14)
pl.savefig('fig%02d.png'%ifig, dpi=dpiset)

#2. 
ifig+=1; pl.figure(ifig)
pl.plot(ob.t, ob.En_imp_eV_time)
pl.xlabel('Time [s]')
pl.ylabel('Deposited electron energy [eV]') 
ms.scix(); pl.grid('on')
pl.suptitle('Var. name: En_imp_eV_time\nElectron energy deposited on the walls in each time-step')
pl.subplots_adjust(top=.82, bottom=.14)
pl.savefig('fig%02d.png'%ifig, dpi=dpiset)

#1. 
ifig+=1; pl.figure(ifig)
pl.plot(ob.t, ob.En_emit_eV_time)
pl.xlabel('Time [s]')
pl.ylabel('Emitted electron energy [eV]') 
ms.scix(); pl.grid('on')
pl.suptitle('Var. name: En_emit_eV_time\nElectron energy emitted by the walls in each time-step')
pl.subplots_adjust(top=.82, bottom=.14)
pl.savefig('fig%02d.png'%ifig, dpi=dpiset)

#3.
ifig+=1; pl.figure(ifig)
pl.plot(ob.t, ob.En_kin_eV_time)
pl.xlabel('Time [s]')
pl.ylabel('$e^-$ kinetic energy [eV]') 
ms.scix(); pl.grid('on')
pl.suptitle('Var. name: En_kin_eV_time\nTotal kinetic energy of the electrons at each time-step')
pl.subplots_adjust(top=.82, bottom=.14)
pl.savefig('fig%02d.png'%ifig, dpi=dpiset)

#7.
ifig+=1; pl.figure(ifig)
pl.plot(ob.t, ob.cen_density)
pl.xlabel('Time [s]')
pl.ylabel('$e^-$ density [$m^{-3}$]')
ms.scix(); pl.grid('on')
pl.suptitle('Var. name: cen_density\nelectron density at the beam position')
pl.subplots_adjust(top=.82, bottom=.14)
pl.savefig('fig%02d.png'%ifig, dpi=dpiset)

#The output file also contains variables represented by matrices with dimension (t_hist, xg_hist) where t_hist[i] represent the time
#right before the i-th passage in the machine and xg_hist[i] represents the position of the i-th slice inside the chamber
#--> dim(t_hist)=#passages, dim(xg_hist)=#slices 



pl.show()






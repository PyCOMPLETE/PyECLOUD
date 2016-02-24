import myloadmat_to_obj as mlo
import pylab as pl
import numpy as np

#load simulation output file 
ob=mlo.myloadmat_to_obj('Pyecltest.mat')

#plot output variables

#VARIABLES VS TIME


#3.Kinetic energy of e-
pl.figure(3)
pl.plot(ob.t, ob.En_kin_eV_time)
pl.xlabel('Time [s]')
pl.ylabel('$e^-$ kinetic energy [eV]') 
pl.title('En_kin_eV_time')
pl.savefig('variable3.png', dpi=300)


#7.Density of electrons inside the chamber
pl.figure(7)
pl.plot(ob.t, ob.cen_density)
pl.xlabel('Time [s]')
pl.ylabel('$e^-$ density [$m^{-3}$]')
pl.title('cen_density') 
pl.savefig('variable7.png', dpi=300)


#8.Beam profile
pl.figure(8)
pl.plot(ob.t, ob.lam_t_array)
pl.xlabel('Time [s]')
pl.ylabel('Beam profile [p/m]') 
pl.title('lam_t_array')
pl.savefig('variable8.png', dpi=300)

#The output file also contains variables represented by matrices with dimension (t_hist, xg_hist) where t_hist[i] represent the time
#corresponding to the i-th passage in the machine and xg_hist[i] represents the position of the i-th slice inside the chamber
#--> dim(t_hist)=#passages, dim(xg_hist)=#slices 



#P |----------------------|
#A |----------------------|
#S |----------------------|
#S |----------------------|
#A |----------------------|
#G |----------------------|
#E |----------------------|
			#SLICE



#Examples
#9.Number of electrons for each passage
pl.figure(9)
pl.plot(np.sum(ob.nel_hist, axis=1)) #axis=1: sum w.r.t. columns 
pl.xlabel('Passage')
pl.ylabel('Number of $e^-$ per passage') 
pl.title('nel_hist')
pl.savefig('variable9.png', dpi=300)


#10.Number of electrons for each slice
pl.figure(10)
pl.plot(ob.xg_hist, np.sum(ob.nel_hist, axis=0)) #axis=0: sum w.r.t. rows
pl.xlabel('Position in the chamber [m]')
pl.ylabel('Number of $e^-$ per slice') 
pl.title('nel_hist')
pl.savefig('variable10.png', dpi=300)

#11.Number of electrons that impact for each passage after scrubbing
pl.figure(11)
pl.plot(np.sum(ob.nel_impact_hist_scrub, axis=1)) 
pl.xlabel('Passage')
pl.ylabel('Number of $e^-$ per passage after scrubbing')
pl.title('nel_impact_hist_scrub') 
pl.savefig('variable11.png', dpi=300)

#12.Average electrons that impact for each slice after scrubbing
pl.figure(12)
pl.plot(ob.xg_hist, np.mean(ob.nel_impact_hist_scrub, axis=0)) 
pl.xlabel('Position in the chamber [m]')
pl.ylabel('Average $e^-$ per slice after scrubbing') 
pl.title('nel_impact_hist_scrub') 
pl.savefig('variable12.png', dpi=300)

#13.Number of electrons that impact for each passage 
pl.figure(13)
pl.plot(np.sum(ob.nel_impact_hist_tot, axis=1)) 
pl.xlabel('Passage')
pl.ylabel('Number of $e^-$ per passage') 
pl.title('nel_impact_hist_tot') 
pl.savefig('variable13.png', dpi=300)

#14.Average electrons that impact for each slice 
pl.figure(14)
pl.plot(ob.xg_hist, np.mean(ob.nel_impact_hist_tot, axis=0)) 
pl.xlabel('Position in the chamber [m]')
pl.ylabel('Average $e^-$ per slice') 
pl.title('nel_impact_hist_tot') 
pl.savefig('variable14.png', dpi=300)

#15.Number of macroparticles for each passage 
pl.figure(15)
pl.plot(ob.N_mp_pass) 
pl.xlabel('Passage')
pl.ylabel('Number of MP per passage') 
pl.title('N_mp_pass') 
pl.savefig('variable15.png', dpi=300)


#16.Number of macroparticles that impact for each passage 
pl.figure(16)
pl.plot(ob.N_mp_impact_pass) 
pl.xlabel('Passage')
pl.ylabel('Number of impacting MP per passage') 
pl.title('N_mp_impact_pass')
pl.savefig('variable16.png', dpi=300)

#17.Number of corrected macroparticles for each passage 
pl.figure(17)
pl.plot(ob.N_mp_corrected_pass) 
pl.xlabel('Passage')
pl.ylabel('Number of corrected MP per passage') 
pl.title('N_mp_corrected_pass')
pl.savefig('variable17.png', dpi=300)


#18.Reference macroparticle size for each passage (this number is dinamically adapted during simulation)
pl.figure(18)
pl.plot(ob.N_mp_ref_pass) 
pl.xlabel('Passage')
pl.ylabel('Reference MP size') 
pl.title('N_mp_ref_pass')
pl.savefig('variable18.png', dpi=300)


#19.Energy of impacting electrons for each passage
pl.figure(19)
pl.plot(np.sum(ob.energ_eV_impact_hist, axis=1))
pl.xlabel('Passage')
pl.ylabel('Energy of impacting electrons for each passage [eV]') 
pl.title('energ_eV_impact_hist')
pl.savefig('variable19.png', dpi=300)

#20. Average energy of impacting electrons for each slice
pl.figure(20)
pl.plot(ob.xg_hist, np.mean(ob.energ_eV_impact_hist, axis=0))
pl.xlabel('Position in the chamber [m]')
pl.ylabel('Average energy of impacting electrons for each slice [eV]') 
pl.title('energ_eV_impact_hist')
pl.savefig('variable20.png', dpi=300)

#21.Electrostatic energy of electrons in the chamber 
pl.figure (21)
pl.plot(ob.t_sc_video, ob.U_sc_eV)
pl.xlabel('Time[s]')
pl.ylabel('Electrostatic energy [eV]') 
pl.title('U_sc_eV')
pl.savefig('variable21.png', dpi=300)


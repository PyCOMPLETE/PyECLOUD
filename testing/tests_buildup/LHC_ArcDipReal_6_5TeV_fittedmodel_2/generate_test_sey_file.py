import sys, os
BIN = os.path.expanduser("../../../../../PyECLOUD/")
sys.path.append(BIN)

costheta_gen = 1.
costheta_test = 0.1#.5
Emax = 332. 
del_max = 1.6
R0 = .7
E0 = 150.



E_start = 0.
E_end = 5000.
DE = 1.

#@profile
def test():

	fname = 'sey_curves.txt'

	import sec_emission_model_ECLOUD as sem_ECLOUD
	import sec_emission_model_FromFile as sem_FromFile
	import numpy as np
	import pylab as pl
	E_vect=np.arange(E_start, E_end, DE)


	delta_vec_ecl, ref_frac_vec_ecl=sem_ECLOUD.yield_fun2(E_vect,costheta_gen+0*E_vect,Emax,del_max,R0,E0)
	
	
	
	pl.close('all')
	pl.figure(1)
	pl.plot(E_vect, delta_vec_ecl)

	data_to_file = np.array([E_vect, ref_frac_vec_ecl*delta_vec_ecl, delta_vec_ecl]).T
	np.savetxt(fname, data_to_file, fmt='%.6e', delimiter=' ', newline='\n', header='', footer='', comments='# ')



	#~ sey_mod=SEY_model_ECLOUD(Emax,del_max,R0)
	sem_ECLOUD_obj = sem_ECLOUD.SEY_model_ECLOUD(Emax=Emax, del_max=del_max, R0=R0 ,E0=E0)
	nel_emit_ec, flag_elast_ec, flag_truesec_ec = sem_ECLOUD_obj.SEY_process(nel_impact=1.+E_vect*0., E_impact_eV=E_vect, costheta_impact=costheta_test+0*E_vect, i_impact=None)

	sem_FromFile_obj = sem_FromFile.SEY_model_FromFile(fname, pyecl_input_folder='./',RThetaDependance=False)
	nel_emit_ff, flag_elast_ff, flag_truesec_ff = sem_FromFile_obj.SEY_process(nel_impact=1.+E_vect*0., E_impact_eV=E_vect, costheta_impact=costheta_test+0*E_vect, i_impact=None)

	
	
	
	pl.plot(E_vect, nel_emit_ec, 'o',color='blue')
	pl.plot(E_vect, nel_emit_ff, 'x',color='red')
	pl.xscale('log')
	N_tests = 10000

	elast_ec =  nel_emit_ff*0.
	elast_ff =  nel_emit_ff*0.

	truesec_ec = nel_emit_ff*0.
	truesec_ff = nel_emit_ff*0.

	for i_test in xrange(N_tests):
		nel_emit_ec, flag_elast_ec, flag_truesec_ec = sem_ECLOUD_obj.SEY_process(nel_impact=1.+E_vect*0., E_impact_eV=E_vect, costheta_impact=costheta_test+0*E_vect, i_impact=None)
		nel_emit_ff, flag_elast_ff, flag_truesec_ff = sem_FromFile_obj.SEY_process(nel_impact=1.+E_vect*0., E_impact_eV=E_vect, costheta_impact=costheta_test+0*E_vect, i_impact=None)
		
		elast_ec += flag_elast_ec
		elast_ff += flag_elast_ff
		
		truesec_ec += flag_truesec_ec
		truesec_ff += flag_truesec_ff

	pl.figure(2)
	pl.plot(E_vect, elast_ec/N_tests)
	pl.plot(E_vect, elast_ff/N_tests)
	pl.xscale('log')
	pl.figure(3)
	pl.plot(E_vect, truesec_ec/N_tests)
	pl.plot(E_vect, truesec_ff/N_tests)
	pl.xscale('log')
	
	pl.show()
	
test()




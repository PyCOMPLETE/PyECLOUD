cimport boris_cython
import numpy as np
cimport numpy as np


cpdef boris_step_multipole(N_sub_steps, Dtt, np.ndarray B_multip, 
						   np.ndarray xn1, np.ndarray yn1, np.ndarray zn1, 
						   np.ndarray vxn1, np.ndarray vyn1, np.ndarray vzn1, 
						   np.ndarray Ex_n, np.ndarray Ey_n, charge, mass):
	
	
	cdef double* B_multip_data = <double*>B_multip.data
	cdef double* xn1_data =  <double*>xn1.data
	cdef double* yn1_data =  <double*>yn1.data
	cdef double* zn1_data =  <double*>zn1.data
	cdef double* vxn1_data =  <double*>vxn1.data
	cdef double* vyn1_data =  <double*>vyn1.data
	cdef double* vzn1_data =  <double*>vzn1.data
	cdef double* Ex_n_data =  <double*>Ex_n.data
	cdef double* Ey_n_data =  <double*>Ey_n.data
	
	boris_c(N_sub_steps, Dtt, B_multip_data, 
		  xn1_data, yn1_data,  zn1_data, 
		  vxn1_data, vyn1_data, vzn1_data,
		  Ex_n_data, Ey_n_data, len(xn1), len(B_multip),
		  charge, mass)


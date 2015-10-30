cdef extern from "boris_c_function.h":
	
	void boris_c(int N_sub_steps, double Dtt, double* B_multip, 
						  double* xn1, double* yn1,  double* zn1, 
						  double* vxn1, double* vyn1, double* vzn1,
						  double* Ex_n, double* Ey_n, int N_mp, int N_multipoles)


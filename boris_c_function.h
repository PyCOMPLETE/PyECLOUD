#ifndef __BORISCFUN
#define __BORISCFUN

void boris_c(int N_sub_steps, double Dtt, 
		double* B_field, double* B_skew,
		double* xn1, double* yn1,  double* zn1, 
		double* vxn1, double* vyn1, double* vzn1,
		double* Ex_n, double* Ey_n, int N_mp, int N_multipoles,
		double charge, double mass);

#endif

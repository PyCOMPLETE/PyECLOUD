#ifndef __BORISCFUN
#define __BORISCFUN


void boris_c(int N_sub_steps, double Dtt, 
		double* k_multip, double* k_skew,  
		double* xn1, double* yn1,  double* zn1, 
		double* vxn1, double* vyn1, double* vzn1,
		double* Ex_n, double* Ey_n, int N_mp, int N_multipoles);


float calc_bx(int order, float param_norm, float param_skew, float x, float y);
float calc_by(int order, float param_norm, float param_skew, float x, float y);


#endif

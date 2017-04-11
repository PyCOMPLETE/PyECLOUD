#ifndef __BORISCFUN
#define __BORISCFUN

#define NONE 		0
#define QUADRUPOLE 	1
#define GENERAL 	2


void boris_c(int N_sub_steps, double Dtt, 
		double* B_multip, double* B_skew,  
		double* xn1, double* yn1,  double* zn1, 
		double* vxn1, double* vyn1, double* vzn1,
		double* Ex_n, double* Ey_n, int N_mp, int N_multipoles);

void b_field_quadrupole(double* B_multip, double xn, double yn, double *Bx, double *By);
void b_field_general(double* B_multip, double* B_skew, double xn, double yn, double *Bx, double *By, int N_multipoles);

int get_b_field_function(double* B_multip, double* B_skew, int N_multipoles, double* Bx, double* By);

#endif

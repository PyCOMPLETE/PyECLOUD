#include "boris_c_function.h"
#include <stdio.h>

const double me = 9.10938291e-31;
const double qe = 1.602176565e-19;

void b_field_drift(double* B_multip, double* B_skew, double xn, double yn, double *Bx, double *By, int N_multipoles){
}

void b_field_dipole(double* B_multip, double* B_skew, double xn, double yn, double *Bx, double *By, int N_multipoles){
	*By = B_multip[0];
}

void b_field_quadrupole(double* B_multip, double* B_skew, double xn, double yn, double *Bx, double *By, int N_multipoles){
	*By = B_multip[1] * xn;
	*Bx = B_multip[1] * yn;
}

void b_field_general(double* B_multip, double* B_skew, double xn, double yn, double *Bx, double *By, int N_multipoles){
	double B_mul_curr;
	double B_skew_curr;
	double rexy_0 ;
	double rexy = 1.;
	double imxy = 0.;

	*By = B_multip[0];
	*Bx = B_skew[0];

	//Order=1 for quadrupoles
	for(int order = 1; order < N_multipoles; order++){
		rexy_0 = rexy;
		rexy = rexy_0*xn - imxy*yn;
		imxy = imxy*xn + rexy_0*yn;
		B_mul_curr = B_multip[order];
		B_skew_curr = B_skew[order];
		*By += (B_mul_curr * rexy - B_skew_curr * imxy);
		*Bx += (B_mul_curr * imxy + B_skew_curr * rexy);
	}
}

b_field_function get_b_field(double* B_multip, double* B_skew, int N_multipoles){
	for(int i=0; i < N_multipoles; i++){
		printf("%.2f\t%.2f\n", B_multip[i], B_skew[i]); 
		if (B_skew[i] != 0.){
			printf("General\n");
			return b_field_general;
		}
	}
	if (N_multipoles == 1){
		if (B_multip[0] == 0.){
			printf("Drift\n");
			return b_field_drift;
		} else {
			printf("Dipole\n");
			return b_field_dipole;
		}
	}
	else if (N_multipoles == 2 && B_multip[0] == 0.){
		printf("Quadrupole\n");
		return b_field_quadrupole;
	}
	else {
		printf("General\n");
		return b_field_general ;
	}
}

void boris_c(int N_sub_steps, double Dtt,
		double* B_multip, double* B_skew,
		double* xn1, double* yn1,  double* zn1,
		double* vxn1, double* vyn1, double* vzn1,
		double* Ex_n, double* Ey_n, int N_mp, int N_multipoles)
{
	int p, isub;
	double Ex_np, Ey_np;
	double Bx_n = 0.;
	double By_n = 0.;

	double tBx, tBy, tBsq;
	double sBx, sBy;
	double vx_prime, vy_prime, vz_prime;
	double vx_min, vy_min, vz_min;
	double vx_plus, vy_plus, vz_plus;
	double vxn1p, vyn1p, vzn1p;
	double xn1p, yn1p, zn1p;
	const double qm = -qe/me; //is an electron

	b_field_function b_field = get_b_field(B_multip, B_skew, N_multipoles);


	for(p=0; p<N_mp; p++)
	{
		Ex_np = Ex_n[p];
		Ey_np = Ey_n[p];

		vxn1p = vxn1[p];
		vyn1p = vyn1[p];
		vzn1p = vzn1[p];

		xn1p = xn1[p];
		yn1p = yn1[p];
		zn1p = zn1[p];

		for (isub=0; isub<N_sub_steps; isub++)
		{
			b_field(B_multip, B_skew, xn1p, yn1p, &Bx_n, &By_n, N_multipoles);

			if (p == 0 && isub == 0){
				printf("%.2f\t%.2f\n", Bx_n, By_n);
			}

			tBx = 0.5*qm*Dtt*Bx_n;
			tBy = 0.5*qm*Dtt*By_n;
			tBsq = tBx*tBx + tBy*tBy; 

			sBx = 2.*tBx/(1.+tBsq);
			sBy = 2.*tBy/(1.+tBsq);

			vx_min = vxn1p + 0.5*qm*Ex_np*Dtt;
			vy_min = vyn1p + 0.5*qm*Ey_np*Dtt;
			vz_min = vzn1p;

			//v_prime = v_min + cross(v_min, tB)
			vx_prime = -vz_min*tBy + vx_min; 
			vy_prime = vz_min*tBx + vy_min;
			vz_prime = vx_min*tBy-vy_min*tBx + vz_min;

			//v_plus = v_min + cross(v_prime, sB)
			vx_plus = -vz_prime*sBy + vx_min;
			vy_plus = vz_prime*sBx + vy_min;
			vz_plus = vx_prime*sBy-vy_prime*sBx + vz_min;

			vxn1p = vx_plus + 0.5*qm*Ex_np*Dtt;
			vyn1p = vy_plus + 0.5*qm*Ey_np*Dtt;
			vzn1p = vz_plus;
			
			xn1p = xn1p + vxn1p * Dtt;
			yn1p = yn1p + vyn1p * Dtt;
			zn1p = zn1p + vzn1p * Dtt;
		}

		xn1[p] = xn1p; 
		yn1[p] = yn1p;
		zn1[p] = zn1p;

		vxn1[p] = vxn1p;
		vyn1[p] = vyn1p;
		vzn1[p] = vzn1p;
	}
}

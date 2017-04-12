#include "boris_c_function.h"

const double me = 9.10938291e-31;
const double qe = 1.602176565e-19;

void boris_c(int N_sub_steps, double Dtt,
		double* B_multip, double* B_skew,
		double* xn1, double* yn1,  double* zn1,
		double* vxn1, double* vyn1, double* vzn1,
		double* Ex_n, double* Ey_n, int N_mp, int N_multipoles)
{
	int p, isub, order;
	double Ex_np, Ey_np;
	double Bx_n, By_n;
	double rexy, imxy, rexy_0;
	double tBx, tBy, tBsq;
	double sBx, sBy;
	double vx_prime, vy_prime, vz_prime;
	double vx_min, vy_min, vz_min;
	double vx_plus, vy_plus, vz_plus;
	double vxn1p, vyn1p, vzn1p;
	double xn1p, yn1p, zn1p;
	const double qm = -qe/me; //is an electron

	/* Sets b-field for drift and dipoles */
	int b_field_type = get_b_field_type(B_multip, B_skew, N_multipoles, &Bx_n, &By_n);

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
			switch( b_field_type ){
				case NONE:
					break;
				case QUADRUPOLE:
					By_n = B_multip[1] * xn1p;
					Bx_n = B_multip[1] * yn1p;
					break;
				case GENERAL:
					rexy = 1.;
					imxy = 0.;
					By_n = B_multip[0];
					Bx_n = B_skew[0];

					for(order = 1; order < N_multipoles; order++){
						/* rexy, imxy correspond to real, imaginary part of (x+iy)^(n-1) */
						rexy_0 = rexy;
						rexy = rexy_0*xn1p - imxy*yn1p;
						imxy = imxy*xn1p + rexy_0*yn1p;

						/*
						 * Bx +iBy = sum[ (k + ik')(x + iy)^(n-1) ]
						 * where k, k' are the strengths and skew strengths of the magnet
						 */
						By_n += (B_multip[order]*rexy - B_skew[order]*imxy);
						Bx_n += (B_multip[order]*imxy + B_skew[order]*rexy);
					}
					break;
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


/*
 * 1. If N_multipoles is 1, the magnetic field is "computed" here and an empty function is returned.
 * 2. Else if it is a non skew quadrupole, a simplified function is returned.
 * 3. Otherwise, the full multipole/skew function is returned.
 */
int get_b_field_type(double* B_multip, double* B_skew, int N_multipoles, double* Bx, double* By){
	//1. Drift or Dipole
	if (N_multipoles == 1){
		*By = B_multip[0];
		*Bx = B_skew[0];
		return NONE;
	}
	//3. Skew quad or skew higher order
	for(int i=0; i < N_multipoles; i++){
		if (B_skew[i] != 0.){
			return GENERAL;
		}
	}
	//2. Simple quad
	if (N_multipoles == 2 && B_multip[0] == 0.){
		return QUADRUPOLE;
	}
	//3. Non skew higher order
	return GENERAL;
}


#include "boris_c_function.h"




void boris_c(int N_sub_steps, double Dtt, double* B_multip, 
						  double* xn1, double* yn1,  double* zn1, 
						  double* vxn1, double* vyn1, double* vzn1,
						  double* Ex_n, double* Ey_n, int N_mp, int N_multipoles)
{
							  
	int p, m, isub;
	double Ex_np, Ey_np;
	double Bx_n, By_n;

	double me, qe, qm;
	double tBx, tBy, tBsq;
	double sBx, sBy;
	double vx_prime, vy_prime, vz_prime;
	double vx_min, vy_min, vz_min;
	double vx_plus, vy_plus, vz_plus;
	double vxn1p, vyn1p, vzn1p;
	double xn1p, yn1p, zn1p;
	double x_to_m, y_to_m;
	double B_mul_curr;
		
	me=9.10938291e-31;
	qe=1.602176565e-19;
	qm=-qe/me; //is an electron

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

			//Compute B field
			
			//Bx_n = 0.;
			//By_n = 0.;
			//x_to_m = 1.;
			//y_to_m = 1.;
			
			//for (m=0;m<N_multipoles;m++)(WRONG! you need the cross terms)
			//{
				//B_mul_curr = B_multip[m];
				//Bx_n = Bx_n + B_mul_curr*y_to_m;
				//By_n = By_n + B_mul_curr*x_to_m;
				//x_to_m = x_to_m*xn1p;
				//y_to_m = y_to_m*yn1p;
			//}
			
			
			//Just dipole and quadrupole for now, to be generalized to multipole
			B_mul_curr = B_multip[0];
			Bx_n = 0.;
			By_n = B_mul_curr;
			
			if (N_multipoles>1)
			{
				B_mul_curr = B_multip[1];
				Bx_n += B_mul_curr*yn1p;
				By_n += B_mul_curr*xn1p;
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





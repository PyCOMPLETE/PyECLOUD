#include "boris_c_function.h"
#include <stdio.h>


int main()
{
	double b_x = 0.;
	double b_y = 0.;

	double xn1p = 0.01;
	double yn1p = 0.00;
	printf("X: %.2f, Y: %.2f\n", xn1p, yn1p);

	double param_norm[] = { 1., 0.7, 1.};
	double param_skew[] = { 0., 0.7, 0.};

	int N_multipoles = 3;

	for(int order=1; order<=N_multipoles; order++)
	{
		double k_mul_curr = param_norm[order-1];
		double k_skew_curr = param_skew[order-1];
		calc_b(order, &k_mul_curr, &k_skew_curr, &xn1p, &yn1p, &b_x, &b_y);

		printf("Bx: %.5f, By: %.5f\n", b_x, b_y);
	}
}

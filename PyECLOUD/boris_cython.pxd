cdef extern from "boris_c_function.h":

    void boris_c(int N_sub_steps, double Dtt, double* B_field, double* B_skew,
                          double* xn1, double* yn1,  double* zn1,
                          double* vxn1, double* vyn1, double* vzn1,
                          double* Ex_n, double* Ey_n,
                          double* Bx_n_custom, double* By_n_custom, double* Bz_n_custom,
		                      int custom_B,
                          int N_mp, int N_multipoles,
                          double charge, double mass)

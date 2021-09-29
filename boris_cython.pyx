cimport boris_cython
import numpy as np
cimport numpy as np


cpdef boris_step_multipole(N_sub_steps, Dtt, np.ndarray B_field, np.ndarray B_skew,
                           np.ndarray xn1, np.ndarray yn1, np.ndarray zn1,
                           np.ndarray vxn1, np.ndarray vyn1, np.ndarray vzn1,
                           np.ndarray Ex_n, np.ndarray Ey_n,
                           np.ndarray Bx_n, np.ndarray By_n, np.ndarray Bz_n,
                           custom_B,
                           charge, mass):


    cdef double* B_field_data = <double*>B_field.data
    cdef double* B_skew_data = <double*>B_skew.data
    cdef double* xn1_data =  <double*>xn1.data
    cdef double* yn1_data =  <double*>yn1.data
    cdef double* zn1_data =  <double*>zn1.data
    cdef double* vxn1_data =  <double*>vxn1.data
    cdef double* vyn1_data =  <double*>vyn1.data
    cdef double* vzn1_data =  <double*>vzn1.data
    cdef double* Ex_n_data =  <double*>Ex_n.data
    cdef double* Ey_n_data =  <double*>Ey_n.data
    cdef double* Bx_n_data =  <double*>Bx_n.data
    cdef double* By_n_data =  <double*>By_n.data
    cdef double* Bz_n_data =  <double*>Bz_n.data

    boris_c(N_sub_steps, Dtt, B_field_data, B_skew_data,
          xn1_data, yn1_data,  zn1_data,
          vxn1_data, vyn1_data, vzn1_data,
          Ex_n_data, Ey_n_data,
          Bx_n_data, By_n_data, Bz_n_data, custom_B,
          len(xn1), len(B_field),
          charge, mass)

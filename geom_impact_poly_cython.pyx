import numpy as np
cimport numpy as np
cimport cython

#~ from cython.parallel import parallel, prange

from libc.math cimport sqrt

ctypedef np.float_t DOUBLE_t
ctypedef np.int_t INT_t

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef impact_point_and_normal(double[::1] x_in, double[::1] y_in, double[::1] z_in,
                              double[::1] x_out, double[::1] y_out, double[::1] z_out,
                              double[::1] Vx, double[::1] Vy,
                              double[::1] Nx, double[::1] Ny, int N_edg, double resc_fac):
    cdef int N_impacts=len(x_in)
    cdef np.ndarray[DOUBLE_t] x_int = np.zeros((N_impacts,),dtype=np.double)
    cdef np.ndarray[DOUBLE_t] y_int = np.zeros((N_impacts,),dtype=np.double)
    cdef np.ndarray[DOUBLE_t] z_int = np.zeros((N_impacts,),dtype=np.double)
    cdef np.ndarray[DOUBLE_t] Nx_int = np.zeros((N_impacts,),dtype=np.double)
    cdef np.ndarray[DOUBLE_t] Ny_int = np.zeros((N_impacts,),dtype=np.double)
    cdef np.ndarray[INT_t] i_found = np.zeros((N_impacts,),dtype=np.int)

    cdef int i_imp, ii, i_found_curr
    cdef double t_min_curr, t_ii, t_border, t_border_min_curr
    cdef int fould_curr
    cdef double x_in_curr, y_in_curr, x_out_curr, y_out_curr, den


    #with nogil, parallel():
    #for i_imp in prange(N_impacts):
    for i_imp in xrange(N_impacts):
        t_min_curr = 1.
        i_found_curr = -1
        fould_curr = False
        x_in_curr = x_in[i_imp]
        y_in_curr = y_in[i_imp]
        x_out_curr = x_out[i_imp]
        y_out_curr = y_out[i_imp]

        for ii in xrange(N_edg):

            den    = ((y_out_curr-y_in_curr)*(Vx[ii+1]-Vx[ii])+(x_in_curr-x_out_curr)*(Vy[ii+1]-Vy[ii]))
            if den == 0.:
                # it is the case when the normal top the segment is perpendicular to the edge
                # the case case overlapping the edge is not possible (this would not allow Pin inside and Pout outside - a point on the edge is condidered outside)
                # the only case left is segment parallel to tue edge => no intersection
                t_border = -2.
            else:
                t_border=((y_out_curr-y_in_curr)*(x_in_curr-Vx[ii])+(x_in_curr-x_out_curr)*(y_in_curr-Vy[ii]))/den


            if t_border>=0. and t_border<=1.:
                t_ii = (Nx[ii]*(Vx[ii]-x_in_curr)+Ny[ii]*(Vy[ii]-y_in_curr)) /(Nx[ii]*(x_out_curr-x_in_curr)+Ny[ii]*(y_out_curr-y_in_curr))
                if t_ii>=0. and t_ii<t_min_curr:
                    t_min_curr=t_ii
                    fould_curr = True
                    i_found_curr = ii


        t_min_curr=resc_fac*t_min_curr
        x_int[i_imp]=t_min_curr*x_out_curr+(1.-t_min_curr)*x_in_curr
        y_int[i_imp]=t_min_curr*y_out_curr+(1.-t_min_curr)*y_in_curr
        z_int[i_imp]=0

        if i_found_curr>=0:
            Nx_int[i_imp] = Nx[i_found_curr]
            Ny_int[i_imp] = Ny[i_found_curr]
            i_found[i_imp] = i_found_curr

    return x_int, y_int, z_int, Nx_int, Ny_int, i_found


ctypedef np.int8_t INT8_t
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef is_outside_convex(np.ndarray[DOUBLE_t] x_mp, np.ndarray[DOUBLE_t] y_mp, np.ndarray[DOUBLE_t] Vx, np.ndarray[DOUBLE_t] Vy, double cx, double cy, int N_edg):

    cdef int N_mp = len(x_mp)
    cdef np.ndarray[INT_t] flag_outside_vec
    cdef int i_mp
    cdef int  ii
    cdef int flag_inside_curr
    cdef double x_curr, y_curr

    #print 'Convex'


    flag_outside_vec = np.zeros((N_mp,),dtype=np.int)

#~     with nogil, parallel():
#~     for i_mp in prange(N_mp):
    for i_mp in xrange(N_mp):
        x_curr = x_mp[i_mp]
        y_curr = y_mp[i_mp]
        flag_inside_curr = (((x_curr/cx)**2 + (y_curr/cy)**2)<=1.)
        #print 1, flag_inside_curr

        if flag_inside_curr==0:
            flag_inside_curr=1
            ii = 0
            while flag_inside_curr==1 and ii<N_edg:
                flag_inside_curr=(((y_curr-Vy[ii])*(Vx[ii+1]-Vx[ii])-(x_curr-Vx[ii])*(Vy[ii+1]-Vy[ii]))>0.)
                ii = ii +1

        flag_outside_vec[i_mp] = not(flag_inside_curr)

    return np.bool_(flag_outside_vec)



@cython.boundscheck(False)
@cython.wraparound(False)
cpdef is_outside_nonconvex(np.ndarray[DOUBLE_t] x_mp, np.ndarray[DOUBLE_t] y_mp, np.ndarray[DOUBLE_t] Vx, np.ndarray[DOUBLE_t] Vy, double cx, double cy, int N_edg):

    cdef int N_mp = len(x_mp)
    cdef np.ndarray[INT_t] flag_outside_vec
    cdef int i_mp
    cdef int  ii, jj
    cdef int flag_inside_curr
    cdef double x_curr, y_curr

    #print 'NON convex'


    flag_outside_vec = np.zeros((N_mp,),dtype=np.int)


    for i_mp in xrange(N_mp):
        x_curr = x_mp[i_mp]
        y_curr = y_mp[i_mp]
        flag_inside_curr = (((x_curr/cx)**2 + (y_curr/cy)**2)<=1.)
        #print 1, flag_inside_curr

        if flag_inside_curr==0:
            ii = 0
            jj = N_edg-1
            while ii < N_edg:
                if ((Vy[ii]>y_curr) != (Vy[jj]>y_curr)):
                    if (x_curr < (Vx[jj]-Vx[ii]) * (y_curr-Vy[ii]) / (Vy[jj]-Vy[ii]) + Vx[ii]) :
                        flag_inside_curr = not(flag_inside_curr)

                jj = ii
                ii += 1


        flag_outside_vec[i_mp] = not(flag_inside_curr)

    return np.bool_(flag_outside_vec)




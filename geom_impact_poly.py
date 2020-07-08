#-Begin-preamble-------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 8.5.0
#
#
#     Main author:          Giovanni IADAROLA
#                           BE-ABP Group
#                           CERN
#                           CH-1211 GENEVA 23
#                           SWITZERLAND
#                           giovanni.iadarola@cern.ch
#
#     Contributors:         Eleonora Belli
#                           Philipp Dijkstal
#                           Lorenzo Giacomel
#                           Lotta Mether
#                           Annalisa Romano
#                           Giovanni Rumolo
#                           Eric Wulff
#
#
#     Copyright  CERN,  Geneva  2011  -  Copyright  and  any   other
#     appropriate  legal  protection  of  this  computer program and
#     associated documentation reserved  in  all  countries  of  the
#     world.
#
#     Organizations collaborating with CERN may receive this program
#     and documentation freely and without charge.
#
#     CERN undertakes no obligation  for  the  maintenance  of  this
#     program,  nor responsibility for its correctness,  and accepts
#     no liability whatsoever resulting from its use.
#
#     Program  and documentation are provided solely for the use  of
#     the organization to which they are distributed.
#
#     This program  may  not  be  copied  or  otherwise  distributed
#     without  permission. This message must be retained on this and
#     any other authorized copies.
#
#     The material cannot be sold. CERN should be  given  credit  in
#     all references.
#
#-End-preamble---------------------------------------------------------


from numpy import squeeze, array, diff, max, sum, sqrt,\
    logical_and, logical_or, ones, zeros, take, arctan2, sin, cos
import scipy.io as sio


class polyg_cham_geom_object:
    def __init__(self, filename_chm, flag_non_unif_sey,
                 flag_verbose_file=False, flag_verbose_stdout=False):

        print('Polygonal chamber - numpy implementation')

        if type(filename_chm) == str:
            dict_chm = sio.loadmat(filename_chm)
        else:
            dict_chm = filename_chm
        Vx = squeeze(dict_chm['Vx'])
        Vy = squeeze(dict_chm['Vy'])
        cx = float(squeeze(dict_chm['x_sem_ellip_insc']))
        cy = float(squeeze(dict_chm['y_sem_ellip_insc']))

        if flag_non_unif_sey == 1:
            self.del_max_segments = squeeze(dict_chm['del_max_segments'])
            self.R0_segments = squeeze(dict_chm['R0_segments'])
            self.Emax_segments = squeeze(dict_chm['Emax_segments'])

        self.N_vert = len(Vx)

        N_edg = len(Vx)

        Vx = list(Vx)
        Vy = list(Vy)

        Vx.append(Vx[0])
        Vy.append(Vy[0])

        Vx = array(Vx)
        Vy = array(Vy)

        Nx = -diff(Vy, 1)
        Ny = diff(Vx, 1)

        norm_N = sqrt(Nx**2 + Ny**2)
        Nx = Nx / norm_N
        Ny = Ny / norm_N

        self.x_aper = max(abs(Vx))
        self.y_aper = max(abs(Vy))
        self.Vx = Vx
        self.Vy = Vy
        self.Nx = Nx
        self.Ny = Ny
        self.N_edg = N_edg
        self.cx = cx
        self.cy = cy

        self.N_mp_impact = 0
        self.N_mp_corrected = 0
        self.chamb_type = 'polyg'

        self.flag_verbose_stdout = flag_verbose_stdout
        self.flag_verbose_file = flag_verbose_file

        if self.flag_verbose_file:
            fbckt = open('bcktr_errors.txt', 'w')
            fbckt.write('kind,x_in,y_in,x_out, y_out\n')
            fbckt.close()

#    def is_outside(self, x_mp, y_mp):
#        N_pts=len(x_mp)
#        flag_inside=array(N_pts*[True])
#        for ii in xrange(self.N_edg):
#            flag_inside[flag_inside]=((y_mp[flag_inside]-self.Vy[ii])*(self.Vx[ii+1]-self.Vx[ii])\
#                                      -(x_mp[flag_inside]-self.Vx[ii])*(self.Vy[ii+1]-self.Vy[ii]))>0
#
#        return ~flag_inside

    def is_outside(self, x_mp, y_mp):

        flag_outside = (((x_mp / self.cx)**2 + (y_mp / self.cy)**2) > 1)
        #flag_outside=array([True]*len(x_mp))

        if flag_outside.any():
            x_mp_chk = x_mp[flag_outside]
            y_mp_chk = y_mp[flag_outside]
            N_pts = len(x_mp_chk)
            flag_inside_chk = array(N_pts * [True])
            for ii in range(self.N_edg):
                flag_inside_chk[flag_inside_chk] = ((y_mp_chk[flag_inside_chk] - self.Vy[ii]) * (self.Vx[ii + 1] - self.Vx[ii])\
                                                    - (x_mp_chk[flag_inside_chk] - self.Vx[ii]) * (self.Vy[ii + 1] - self.Vy[ii])) > 0
            flag_outside[flag_outside] = ~flag_inside_chk

        return flag_outside

    def impact_point_and_normal(self, x_in, y_in, z_in, x_out, y_out, z_out, resc_fac=0.99, flag_robust=True):

        N_impacts = len(x_in)
        t_min = ones(N_impacts) + 1.
        i_found = array(N_impacts * [-1])
        mask_found = array(N_impacts * [False])

        self.N_mp_impact = self.N_mp_impact + N_impacts

        for ii in range(self.N_edg):
            t_curr = (self.Nx[ii] * (self.Vx[ii] - x_in) + self.Ny[ii] * (self.Vy[ii] - y_in)) / \
                (self.Nx[ii] * (x_out - x_in) + self.Ny[ii] * (y_out - y_in))

            mask_updatemin = logical_and(t_curr >= 0, t_curr < t_min)
            t_min[mask_updatemin] = t_curr[mask_updatemin]
            mask_found = logical_or(mask_updatemin, mask_found)
            i_found[mask_updatemin] = ii

            mask_same_min = (t_curr == t_min)
            if mask_same_min.any():
                mask_upd_i_found = array(N_impacts * [False])
                t_border = ((y_out[mask_same_min] - y_in[mask_same_min]) * (x_in[mask_same_min] - self.Vx[ii]) + (x_in[mask_same_min] - x_out[mask_same_min]) * (y_in[mask_same_min] - self.Vy[ii])) / \
                    ((y_out[mask_same_min] - y_in[mask_same_min]) * (self.Vx[ii + 1] - self.Vx[ii]) + (x_in[mask_same_min] - x_out[mask_same_min]) * (self.Vy[ii + 1] - self.Vy[ii]))
                mask_upd_i_found[mask_same_min] = logical_and(t_border >= 0., t_border <= 1.)
                i_found[mask_upd_i_found] = ii

        t_min = resc_fac * t_min
        x_int = t_min * x_out + (1. - t_min) * x_in
        y_int = t_min * y_out + (1. - t_min) * y_in
        z_int = 0 * t_min

        Nx_int = zeros(N_impacts)
        Ny_int = zeros(N_impacts)
        Nx_int[mask_found] = take(self.Nx, i_found[mask_found])
        Ny_int[mask_found] = take(self.Ny, i_found[mask_found])

        if sum(mask_found) < N_impacts:
            mask_not_found = ~mask_found

            x_int[mask_not_found] = x_in[mask_not_found]
            y_int[mask_not_found] = y_in[mask_not_found]

            #compute some kind of normal ....
            par_cross = arctan2(self.cx * y_in[mask_not_found], self.cy * x_int[mask_not_found])

            Dx = -self.cx * sin(par_cross)
            Dy = self.cy * cos(par_cross)

            Nx_corr = -Dy
            Ny_corr = Dx

            neg_flag = ((Nx_corr * x_int[mask_not_found] + Ny_corr * y_int[mask_not_found]) > 0)

            Nx_corr[neg_flag] = -Nx_corr[neg_flag]
            Ny_corr[neg_flag] = -Ny_corr[neg_flag]

            Nx_int[mask_not_found] = Nx_corr
            Ny_int[mask_not_found] = Ny_corr

            x_in_error = x_in[mask_not_found]
            y_in_error = y_in[mask_not_found]
            x_out_error = x_out[mask_not_found]
            y_out_error = y_out[mask_not_found]
            N_errors = len(x_in_error)
            self.N_mp_corrected = self.N_mp_corrected + N_errors

            if self.flag_verbose_stdout:
                print("""Reporting backtrack error of kind 1: no impact found""")
                print("""x_in, y_in, x_out, y_out""")
                for i_err in range(N_errors):
                    lcurr = '%.10e,%.10e,%.10e,%.10e'%(x_in_error[i_err], y_in_error[i_err], x_out_error[i_err], y_out_error[i_err])
                    print(lcurr)
                print("""End reporting backtrack error of kind 1""")

            if self.flag_verbose_file:
                with open('bcktr_errors.txt', 'a') as fbckt:
                    for i_err in range(N_errors):
                        lcurr = '%.10e,%.10e,%.10e,%.10e'%(x_in_error[i_err], y_in_error[i_err], x_out_error[i_err], y_out_error[i_err])
                        fbckt.write('1,' + lcurr + '\n')

        if flag_robust:
            flag_impact = self.is_outside(x_int, y_int)
            if flag_impact.any():
                self.N_mp_corrected = self.N_mp_corrected + sum(flag_impact)
                x_int[flag_impact] = x_in[flag_impact]
                y_int[flag_impact] = y_in[flag_impact]
                x_in_error = x_in[flag_impact]
                y_in_error = y_in[flag_impact]
                x_out_error = x_out[flag_impact]
                y_out_error = y_out[flag_impact]
                N_errors = len(x_in_error)

                if self.flag_verbose_stdout:
                    print("""Reporting backtrack error of kind 2: outside after backtracking""")
                    print("""x_in, y_in, x_out, y_out""")
                    for i_err in range(N_errors):
                        lcurr = '%.10e,%.10e,%.10e,%.10e'%(x_in_error[i_err], y_in_error[i_err], x_out_error[i_err], y_out_error[i_err])
                        print(lcurr)
                    print("""End reporting backtrack error of kind 2""")

                if self.flag_verbose_file:
                    with open('bcktr_errors.txt', 'a') as fbckt:
                        for i_err in range(N_errors):
                            lcurr = '%.10e,%.10e,%.10e,%.10e'%(x_in_error[i_err], y_in_error[i_err], x_out_error[i_err], y_out_error[i_err])
                            fbckt.write('2,' + lcurr + '\n')

            flag_impact = self.is_outside(x_int, y_int)
            if sum(flag_impact) > 0:
                #~ import pylab as pl
                #~ pl.close('all')
                #~ pl.plot(self.Vx, self.Vy)
                #~ pl.plot(x_in, y_in,'.b')
                #~ pl.plot(x_out, y_out,'.k')
                #~ pl.plot(x_int, y_int,'.g')
                #~ pl.plot(x_int[flag_impact], y_int[flag_impact],'.r')
                #~ pl.show()
                if self.flag_verbose_stdout:
                    print("""Reporting backtrack error of kind 3: outside after correction""")
                    print("""x_in, y_in, x_out, y_out""")
                x_in_error = x_in[flag_impact]
                y_in_error = y_in[flag_impact]
                x_out_error = x_out[flag_impact]
                y_out_error = y_out[flag_impact]
                N_errors = len(x_in_error)

                if self.flag_verbose_stdout:
                    print("""Reporting backtrack error of kind 3: outside after correction""")
                    print("""x_in, y_in, x_out, y_out""")
                    for i_err in range(N_errors):
                        lcurr = '%.10e,%.10e,%.10e,%.10e'%(x_in_error[i_err], y_in_error[i_err], x_out_error[i_err], y_out_error[i_err])
                        print(lcurr)
                    print("""End reporting backtrack error of kind 3""")

                if self.flag_verbose_file:
                    with open('bcktr_errors.txt', 'a') as fbckt:
                        for i_err in range(N_errors):
                            lcurr = '%.10e,%.10e,%.10e,%.10e'%(x_in_error[i_err], y_in_error[i_err], x_out_error[i_err], y_out_error[i_err])
                            fbckt.write('3,' + lcurr + '\n')

                raise ValueError('Outside after backtracking!!!!')

        return x_int, y_int, z_int, Nx_int, Ny_int, i_found




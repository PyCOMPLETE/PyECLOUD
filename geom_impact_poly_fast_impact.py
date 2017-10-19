#----------------------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 6.5.1
#
#
#     Author and contact:   Giovanni IADAROLA
#                           BE-ABP Group
#                           CERN
#                           CH-1211 GENEVA 23
#                           SWITZERLAND
#                           giovanni.iadarola@cern.ch
#
#                contact:   Giovanni RUMOLO
#                           BE-ABP Group
#                           CERN
#                           CH-1211 GENEVA 23
#                           SWITZERLAND
#                           giovanni.rumolo@cern.ch
#
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
#----------------------------------------------------------------------


from __future__ import division, print_function
from numpy import squeeze, array,diff, max, sum, sqrt,\
                  arctan2, sin, cos
import scipy.io as sio
import numpy as np
import numpy.random as random

import geom_impact_poly_cython as gipc

class PyECLOUD_ChamberException(ValueError):
    pass

class polyg_cham_geom_object:
    def __init__(self, filename_chm, flag_non_unif_sey, flag_verbose_file=False, flag_verbose_stdout=False,
                 flag_assume_convex=True, flag_counter_clockwise_chamb=True, distance_new_phem=1e-5):


        print('Polygonal chamber - cython implementation')

        if type(filename_chm)==str:
            dict_chm=sio.loadmat(filename_chm)
        else:
            dict_chm=filename_chm
        self.dict_chm = dict_chm

        Vx=squeeze(dict_chm['Vx'])
        Vy=squeeze(dict_chm['Vy'])
        cx=float(squeeze(dict_chm['x_sem_ellip_insc']))
        cy=float(squeeze(dict_chm['y_sem_ellip_insc']))

        if flag_non_unif_sey==1:
            self.del_max_segments = squeeze(dict_chm['del_max_segments'])
            self.R0_segments = squeeze(dict_chm['R0_segments'])
            self.Emax_segments = squeeze(dict_chm['Emax_segments'])

        # For photoemission from segment
        if 'phem_cdf' in dict_chm:
            self.phem_cdf = squeeze(dict_chm['phem_cdf'])
            self.distance_new_phem = distance_new_phem

            # Make sure phem_cdf has correct shape
            if self.phem_cdf[-1] != 1:
                raise PyECLOUD_ChamberException('phem_cdf of chamb_dict does not end with 1.')
            if np.any(np.diff(self.phem_cdf) < 0):
                raise PyECLOUD_ChamberException('phem_cdf of chamb_dict is not monotonically increasing.')

            # Needed to calculate histograms and positions later
            _Vx = np.append(Vx, [Vx[0]])
            _Vy = np.append(Vy, [Vy[0]])
            self.seg_diff_x = seg_diff_x = _Vx[1:] - _Vx[:-1]
            self.seg_diff_y = seg_diff_y = _Vy[1:] - _Vy[:-1]
            self.cdf_bins = np.concatenate(([0], self.phem_cdf))

            _len_segments = np.sqrt(seg_diff_x**2  + seg_diff_y**2)
            self.normal_vect_x = seg_diff_y / _len_segments
            self.normal_vect_y = -seg_diff_x / _len_segments

            # Chamber corners are defined clockwise or counter-clockwise.
            # This has an effect of the direction of the normal vector
            if not flag_counter_clockwise_chamb:
                self.normal_vect_x *= -1
                self.normal_vect_y *= -1

            # MP positions have to be inside chamber -> add normal_vector*distance_new_phem to position
            self.phem_x0 = Vx + self.normal_vect_x*self.distance_new_phem
            self.phem_y0 = Vy + self.normal_vect_y*self.distance_new_phem


        self.N_vert=len(Vx)

        N_edg=len(Vx)

        Vx=list(Vx)
        Vy=list(Vy)

        Vx.append(Vx[0])
        Vy.append(Vy[0])

        Vx=array(Vx)
        Vy=array(Vy)

        Nx=-diff(Vy,1)
        Ny=diff(Vx,1)

        norm_N=sqrt(Nx**2+Ny**2)
        Nx=Nx/norm_N
        Ny=Ny/norm_N

        self.x_aper = max(np.abs(Vx))
        self.y_aper = max(np.abs(Vy))
        self.Vx=Vx
        self.Vy=Vy
        self.Nx=Nx
        self.Ny=Ny
        self.N_edg=N_edg
        self.cx=cx
        self.cy=cy

        self.N_mp_impact=0
        self.N_mp_corrected=0
        self.chamb_type='polyg'

        self.flag_verbose_stdout = flag_verbose_stdout
        self.flag_verbose_file = flag_verbose_file

        self.flag_assume_convex = flag_assume_convex

        if self.flag_verbose_file:
            fbckt=open('bcktr_errors.txt','w')
            fbckt.write('kind,x_in,y_in,x_out, y_out\n')
            fbckt.close()

        if self.flag_assume_convex:
            if not(self.is_convex()):
                raise PyECLOUD_ChamberException(
                    'The polygon looks not convex!!!!\nIn this case you can use the general algorithm (probably slower) by setting:\nflag_assume_convex = False')
            self.cythonisoutside = gipc.is_outside_convex
            print('Assuming convex polygon')
        else:
            self.cythonisoutside = gipc.is_outside_nonconvex
            print('No assumption on the convexity of the polygon')


    def is_outside(self, x_mp, y_mp):
        #~ return gipc.is_outside(x_mp, y_mp, self.Vx, self.Vy, self.cx, self.cy, self.N_edg)
        return self.cythonisoutside(x_mp, y_mp, self.Vx, self.Vy, self.cx, self.cy, self.N_edg)


    #@profile
    def impact_point_and_normal(self, x_in, y_in, z_in, x_out, y_out, z_out, resc_fac=0.99, flag_robust=True):

        N_impacts=len(x_in)
        self.N_mp_impact=self.N_mp_impact+N_impacts

        x_int,y_int,z_int,Nx_int,Ny_int, i_found = gipc.impact_point_and_normal(x_in, y_in, z_in, x_out, y_out, z_out,
                              self.Vx,  self.Vy,self.Nx,  self.Ny,  self.N_edg, resc_fac)

        mask_found = i_found>=0

        if sum(mask_found)<N_impacts:
            mask_not_found = ~mask_found

            x_int[mask_not_found] = x_in[mask_not_found]
            y_int[mask_not_found] = y_in[mask_not_found]

            #compute some kind of normal ....
            par_cross=arctan2(self.cx*y_in[mask_not_found],self.cy*x_int[mask_not_found])

            Dx=-self.cx*sin(par_cross)
            Dy=self.cy*cos(par_cross)

            Nx_corr=-Dy
            Ny_corr=Dx

            neg_flag=((Nx_corr*x_int[mask_not_found]+Ny_corr*y_int[mask_not_found])>0)

            Nx_corr[neg_flag]=-Nx_corr[neg_flag]
            Ny_corr[neg_flag]=-Ny_corr[neg_flag]

            Nx_int[mask_not_found]=Nx_corr
            Ny_int[mask_not_found]=Ny_corr


            x_in_error = x_in[mask_not_found]
            y_in_error = y_in[mask_not_found]
            x_out_error = x_out[mask_not_found]
            y_out_error = y_out[mask_not_found]
            N_errors = len(x_in_error)
            self.N_mp_corrected = self.N_mp_corrected + N_errors

            if self.flag_verbose_stdout:
                print('Reporting backtrack error of kind 1: no impact found')
                print('x_in, y_in, x_out, y_out')
                for i_err in xrange(N_errors):
                    lcurr = '%.10e,%.10e,%.10e,%.10e' % (x_in_error[i_err], y_in_error[i_err], x_out_error[i_err], y_out_error[i_err])
                    print(lcurr)
                print('End reporting backtrack error of kind 1')

            if self.flag_verbose_file:
                with open('bcktr_errors.txt','a') as fbckt:
                    for i_err in xrange(N_errors):
                        lcurr = '%.10e,%.10e,%.10e,%.10e' % (x_in_error[i_err], y_in_error[i_err], x_out_error[i_err], y_out_error[i_err])
                        fbckt.write('1,'+lcurr+'\n')




        if flag_robust:
            flag_impact=self.is_outside(x_int, y_int)
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
                    print('Reporting backtrack error of kind 2: outside after backtracking')
                    print('x_in, y_in, x_out, y_out')
                    for i_err in xrange(N_errors):
                        lcurr = '%.10e,%.10e,%.10e,%.10e' % (x_in_error[i_err], y_in_error[i_err], x_out_error[i_err], y_out_error[i_err])
                        print(lcurr)
                    print('End reporting backtrack error of kind 2')

                if self.flag_verbose_file:
                    with open('bcktr_errors.txt','a') as fbckt:
                        for i_err in xrange(N_errors):
                            lcurr = '%.10e,%.10e,%.10e,%.10e' % (x_in_error[i_err], y_in_error[i_err], x_out_error[i_err], y_out_error[i_err])
                            fbckt.write('2,'+lcurr+'\n')

            flag_impact=self.is_outside(x_int, y_int)
            if sum(flag_impact)>0:
                #~ import pylab as pl
                #~ pl.close('all')
                #~ pl.plot(self.Vx, self.Vy)
                #~ pl.plot(x_in, y_in,'.b')
                #~ pl.plot(x_out, y_out,'.k')
                #~ pl.plot(x_int, y_int,'.g')
                #~ pl.plot(x_int[flag_impact], y_int[flag_impact],'.r')
                #~ pl.show()
                if self.flag_verbose_stdout:
                    print('Reporting backtrack error of kind 3: outside after correction')
                    print('x_in, y_in, x_out, y_out')
                x_in_error = x_in[flag_impact]
                y_in_error = y_in[flag_impact]
                x_out_error = x_out[flag_impact]
                y_out_error = y_out[flag_impact]
                N_errors = len(x_in_error)

                if self.flag_verbose_stdout:
                    print('Reporting backtrack error of kind 3: outside after correction')
                    print('x_in, y_in, x_out, y_out')
                    for i_err in xrange(N_errors):
                        lcurr = '%.10e,%.10e,%.10e,%.10e' % (x_in_error[i_err], y_in_error[i_err], x_out_error[i_err], y_out_error[i_err])
                        print(lcurr)
                    print('End reporting backtrack error of kind 3')

                if self.flag_verbose_file:
                    with open('bcktr_errors.txt','a') as fbckt:
                        for i_err in xrange(N_errors):
                            lcurr = '%.10e,%.10e,%.10e,%.10e' % (x_in_error[i_err], y_in_error[i_err], x_out_error[i_err], y_out_error[i_err])
                            fbckt.write('3,'+lcurr+'\n')


                raise PyECLOUD_ChamberException('Outside after backtracking!!!!')


        return x_int, y_int, z_int, Nx_int, Ny_int, i_found

    def is_convex(self):
            # From:
            # http://csharphelper.com/blog/2014/07/determine-whether-a-polygon-is-convex-in-c/

            # For each set of three adjacent points A, B, C,
            # find the cross product AB x BC. If the sign of
            # all the cross products is the same, the angles
            # are all positive or negative (depending on the
            # order in which we visit them) so the polygon
            # is convex.
            got_negative = False
            got_positive = False
            num_points = self.N_edg
            for A in xrange(num_points):

                B = np.mod((A + 1),num_points)
                C = np.mod((B + 1), num_points)

                BAx = self.Vx[A] - self.Vx[B]
                BAy = self.Vy[A] - self.Vy[B]
                BCx = self.Vx[C] - self.Vx[B]
                BCy = self.Vy[C] - self.Vy[B]

                cross_product = (BAx * BCy - BAy * BCx)

                if (cross_product < 0):
                    got_negative = True
                elif (cross_product > 0):
                    got_positive = True
                if (got_negative and got_positive):
                    return False


            # If we got this far, the polygon is convex.
            return True

    def get_photoelectron_positions(self, N_mp_gen):
        """
        input: N_mp_gen - Number of MPs for which positions and normal vectors are to be calculated
        The cdf that is part of the chamber definition is used to generate the output
        output: positions and normal vectors for every generated macroparticle
        """
        # Only meaningful for photoemission from segment
        x_new_mp = np.empty(N_mp_gen)
        y_new_mp = np.empty(N_mp_gen)
        norm_x_new_mp = np.empty(N_mp_gen)
        norm_y_new_mp = np.empty(N_mp_gen)

        # Divide new MPs into segments according to the cdf
        N_mp_segment, _ = np.histogram(random.rand(N_mp_gen), self.cdf_bins)

        # The MPs in each segment are distributed evenly along the segment
        N_mp_curr = 0
        for i_seg, N_mp_seg in enumerate(N_mp_segment):
            if N_mp_seg != 0:
                N_mp_after = N_mp_curr + N_mp_seg
                rr = random.rand(N_mp_seg)
                x_new_mp[N_mp_curr:N_mp_after] = self.phem_x0[i_seg] + rr*self.seg_diff_x[i_seg]
                y_new_mp[N_mp_curr:N_mp_after] = self.phem_y0[i_seg] + rr*self.seg_diff_y[i_seg]
                norm_x_new_mp[N_mp_curr:N_mp_after] = self.normal_vect_x[i_seg]
                norm_y_new_mp[N_mp_curr:N_mp_after] = self.normal_vect_y[i_seg]
                N_mp_curr = N_mp_after
        assert N_mp_curr == N_mp_gen

        # Some MPs may be outside of the chamber -> recursion until all are inside
        phem_outside = self.is_outside(x_new_mp, y_new_mp)
        n_phem_outside = np.sum(phem_outside)
        if n_phem_outside != 0:
            x_new_mp[phem_outside], y_new_mp[phem_outside], norm_x_new_mp[phem_outside], norm_y_new_mp[phem_outside] =\
                self.get_photoelectron_positions(n_phem_outside)

        return x_new_mp, y_new_mp, norm_x_new_mp, norm_y_new_mp


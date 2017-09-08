#----------------------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 6.3.1
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


from numpy import sqrt, sin, cos, squeeze, sum
import scipy.io as sio
import int_field_for as iff

me=9.10938291e-31;
qe=1.602176565e-19;
qm=qe/me; #formulas are already written for a negative charge


class pusher_strong_B_generalized():
    def __init__(self, Dt, B0x, B0y, \
                 B_map_file, fact_Bmap, B_zero_thrhld):

        print "Tracker: Generalized strong B"

        self.Dt = Dt
        self.B0x = B0x
        self.B0y = B0y
        self.B_zero_thrhld = B_zero_thrhld
        if B_map_file is None:
            self.flag_B_map = False
            self.analyt_quad_grad1 = False
        elif B_map_file is 'analytic_qaudrupole_unit_grad':
            print "B map analytic quadrupole"
            self.flag_B_map = False
            self.analyt_quad_grad1 = True
            self.fact_Bmap = fact_Bmap
        else:
            self.flag_B_map = True
            self.analyt_quad_grad1 = False
            print 'Loading B map'
            dict_Bmap=sio.loadmat(B_map_file)

            self.Bmap_x = fact_Bmap*squeeze(dict_Bmap['Bx'].real)
            self.Bmap_y = fact_Bmap*squeeze(dict_Bmap['By'].real)
            self.xx=squeeze(dict_Bmap['xx'].T)
            self.yy=squeeze(dict_Bmap['yy'].T)

            self.xmin=min(self.xx);
            self.ymin=min(self.yy);
            self.dx=self.xx[1]-self.xx[0];
            self.dy=self.yy[1]-self.yy[0];

            #            ####Debug
#            import pylab as pl
#            pl.figure(1)
#            pl.imshow(self.Bmap_x.T, cmap=None, norm=None, aspect=None, interpolation=None,
#                      alpha=None, vmin=None, vmax=None, origin=None, extent=None)
#            pl.figure(2)
#            pl.imshow(self.Bmap_y.T, cmap=None, norm=None, aspect=None, interpolation=None,
#                      alpha=None, vmin=None, vmax=None, origin=None, extent=None)
#
#            pl.show()
#
#            raise ValueError
#
#            ####



#    def step(self, xn, yn, zn, vxn, vyn, vzn,Ex_n,Ey_n):
    def step(self, MP_e, Ex_n,Ey_n):

        if MP_e.N_mp>0:
            xn = MP_e.x_mp[0:MP_e.N_mp]
            yn = MP_e.y_mp[0:MP_e.N_mp]
            zn = MP_e.z_mp[0:MP_e.N_mp]
            vxn = MP_e.vx_mp[0:MP_e.N_mp]
            vyn = MP_e.vy_mp[0:MP_e.N_mp]
            vzn = MP_e.vz_mp[0:MP_e.N_mp]


            if self.flag_B_map:
                Bx_n,By_n = iff.int_field(xn,yn,self.xmin,self.ymin,\
                                      self.dx,self.dy,self.Bmap_x,self.Bmap_y)
                # the rescaling factor has already been applied to the map
            elif self.analyt_quad_grad1:
                # the rescaling factor has to be applied here
                Bx_n = self.fact_Bmap*yn.copy()
                By_n = self.fact_Bmap*xn.copy()
            else:
                Bx_n = 0*xn
                By_n = 0*xn



            Bx_n = Bx_n + self.B0x
            By_n = By_n + self.B0y

            #Bx_n = 0*xn + self.B0x
            #By_n = 0*xn + self.B0y

            B_mod=sqrt(Bx_n*Bx_n+By_n*By_n)

            flag_B_zero=B_mod<self.B_zero_thrhld

            N_B_zero=sum(flag_B_zero)

            #introduce fake B for zero field points
            if N_B_zero>0:
                B_mod[flag_B_zero]=1.

            omegac=B_mod*qm;
            sinwcDt = sin(omegac*self.Dt)
            coswcDt = cos(omegac*self.Dt)

            vB_Bmod = (vxn*By_n - vyn*Bx_n)/B_mod
            EB_Bmodsq = (Ex_n*By_n - Ey_n*Bx_n)/(B_mod*B_mod)

            Dz =  vzn*sinwcDt/omegac + vB_Bmod*coswcDt/omegac + \
                EB_Bmodsq*(self.Dt - sinwcDt/omegac) - vB_Bmod/omegac

            Dz_tilde = (vzn*(1-coswcDt)+vB_Bmod*sinwcDt)/(omegac*omegac) + \
                       EB_Bmodsq*(self.Dt*self.Dt/2.+(coswcDt-1.)/(omegac*omegac)) \
                       - vB_Bmod*self.Dt/omegac


            vzn1 = vzn*coswcDt - vB_Bmod*sinwcDt + EB_Bmodsq * (1.-coswcDt)
            vxn1 = vxn - qm*(Ex_n*self.Dt - By_n*Dz)
            vyn1 = vyn - qm*(Ey_n*self.Dt + Bx_n*Dz )

            xn1 = xn + vxn*self.Dt - qm * (Ex_n*self.Dt*self.Dt/2. - By_n*Dz_tilde)
            yn1 = yn + vyn*self.Dt - qm * (Ey_n*self.Dt*self.Dt/2. + Bx_n*Dz_tilde)
            zn1 = zn + Dz


            # correcting zero field points
            if N_B_zero>0:

                xn1[flag_B_zero]=xn[flag_B_zero]+vxn[flag_B_zero]*self.Dt\
                    -0.5*qm*Ex_n[flag_B_zero]*self.Dt*self.Dt;
                yn1[flag_B_zero]=yn[flag_B_zero]+vyn[flag_B_zero]*self.Dt\
                    -0.5*qm*Ey_n[flag_B_zero]*self.Dt*self.Dt;
                zn1[flag_B_zero]=zn[flag_B_zero]+vzn[flag_B_zero]*self.Dt

                vxn1[flag_B_zero]=vxn[flag_B_zero]-qm*Ex_n[flag_B_zero]*self.Dt;
                vyn1[flag_B_zero]=vyn[flag_B_zero]-qm*Ey_n[flag_B_zero]*self.Dt;
                vzn1[flag_B_zero]=vzn[flag_B_zero]

            MP_e.x_mp[0:MP_e.N_mp] = xn1
            MP_e.y_mp[0:MP_e.N_mp] = yn1
            MP_e.z_mp[0:MP_e.N_mp] = zn1
            MP_e.vx_mp[0:MP_e.N_mp] = vxn1
            MP_e.vy_mp[0:MP_e.N_mp] = vyn1
            MP_e.vz_mp[0:MP_e.N_mp]  = vzn1



        return MP_e



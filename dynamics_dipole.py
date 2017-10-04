#----------------------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 6.4.1
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

import sincc_cosincc_cubsincc as cs


me=9.10938291e-31;
qe=1.602176565e-19;
qm=qe/me;


class pusher_dipole_magnet():
    def __init__(self, Dt, B):


        print "Tracker: Dipole - strong B"
        omegac=B*qm;
        self.sincc_v=cs.sincc(omegac*Dt);
        self.cosincc_v=cs.cosincc(omegac*Dt);
        self.cubsincc_v=cs.cubsincc(omegac*Dt);
        self.sin_v=cs.sin(omegac*Dt);
        self.cos_v=cs.cos(omegac*Dt);
        self.Dt=Dt
        self.omegac=omegac

#    def step(self, xn, yn, zn, vxn, vyn, vzn,Ex_n,Ey_n):
    def step(self, MP_e, Ex_n,Ey_n):

        if MP_e.N_mp>0:
            xn = MP_e.x_mp[0:MP_e.N_mp]
            yn = MP_e.y_mp[0:MP_e.N_mp]
            zn = MP_e.z_mp[0:MP_e.N_mp]
            vxn = MP_e.vx_mp[0:MP_e.N_mp]
            vyn = MP_e.vy_mp[0:MP_e.N_mp]
            vzn = MP_e.vz_mp[0:MP_e.N_mp]


            xn1=xn+vxn*self.sincc_v*self.Dt+vzn*self.cosincc_v*self.omegac*self.Dt*self.Dt\
                -qm*Ex_n*self.cosincc_v*self.Dt*self.Dt;
            yn1=yn+vyn*self.Dt-0.5*qm*Ey_n*self.Dt*self.Dt;
            zn1=zn+vzn*self.sincc_v*self.Dt-vxn*self.cosincc_v*self.omegac*self.Dt*self.Dt\
                +qm*Ex_n*self.cubsincc_v*self.omegac*self.Dt*self.Dt*self.Dt;

            vxn1=vxn*self.cos_v+vzn*self.sin_v-qm*Ex_n*self.sincc_v*self.Dt;
            vyn1=vyn-qm*Ey_n*self.Dt;
            vzn1=vzn*self.cos_v-vxn*self.sin_v+qm*Ex_n*self.cosincc_v*self.omegac*self.Dt*self.Dt;


            MP_e.x_mp[0:MP_e.N_mp] = xn1
            MP_e.y_mp[0:MP_e.N_mp] = yn1
            MP_e.z_mp[0:MP_e.N_mp] = zn1
            MP_e.vx_mp[0:MP_e.N_mp] = vxn1
            MP_e.vy_mp[0:MP_e.N_mp] = vyn1
            MP_e.vz_mp[0:MP_e.N_mp]  = vzn1



        return MP_e




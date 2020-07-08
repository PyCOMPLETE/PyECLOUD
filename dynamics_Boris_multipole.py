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

import math
import numpy as np
from .boris_cython import boris_step_multipole


class pusher_Boris_multipole():

    def __init__(self, Dt, N_sub_steps=1, B_multip=None, B_skew=None,
        B0x=None, B0y=None, B0z=None):

        self.N_sub_steps = N_sub_steps
        self.Dt = Dt

        if Dt is None or N_sub_steps is None:
            self.Dtt = None
        else:
            self.Dtt = Dt / float(N_sub_steps)

        if B_multip is None or len(B_multip) == 0:
            B_multip = np.array([0.], dtype=float)

        # B_multip are derivatives of B_field
        # B_field are field strengths at x=1 m, y=0
        factorial = np.array([math.factorial(ii) for ii in range(len(B_multip))], dtype=float)
        self.B_field = np.array(B_multip, dtype=float) / factorial
        if B_skew is None:
            self.B_field_skew = np.zeros_like(self.B_field, dtype=float)
        else:
            self.B_field_skew = np.array(B_skew, dtype=float) / factorial

        self.B0x = B0x
        self.B0y = B0y
        self.B0z = B0z
        print("Tracker: Boris multipole")

        print("N_subst_init=%d" % self.N_sub_steps)

    #@profile
    def step(self, MP_e, Ex_n, Ey_n, Ez_n=0., Bx_n=None, By_n=None, Bz_n=None):
        MP_e = self.stepcustomDt(MP_e, Ex_n, Ey_n, Ez_n,
            Bx_n, By_n, Bz_n,
            Dt_substep=self.Dtt, N_sub_steps=self.N_sub_steps)
        return MP_e

    def stepcustomDt(self, MP_e, Ex_n, Ey_n, Ez_n=0.,
        Bx_n=None, By_n=None, Bz_n=None,
        Dt_substep=None, N_sub_steps=None):

        custom_B = 0
        Bx_arr = np.zeros(MP_e.N_mp)
        By_arr = np.zeros(MP_e.N_mp)
        Bz_arr = np.zeros(MP_e.N_mp)

        if self.B0x is not None:
            Bx_arr += self.B0x
            custom_B = 1
        if self.B0y is not None:
            By_arr += self.B0y
            custom_B = 1
        if self.B0z is not None:
            Bz_arr += self.B0z
            custom_B = 1

        if Bx_n is not None:
            Bx_arr += Bx_n
            custom_B = 1
        if By_n is not None:
            By_arr += By_n
            custom_B = 1
        if Bz_n is not None:
            Bz_arr += Bz_n
            custom_B = 1

        if MP_e.N_mp > 0:

            xn1 = MP_e.x_mp[0:MP_e.N_mp]
            yn1 = MP_e.y_mp[0:MP_e.N_mp]
            zn1 = MP_e.z_mp[0:MP_e.N_mp]
            vxn1 = MP_e.vx_mp[0:MP_e.N_mp]
            vyn1 = MP_e.vy_mp[0:MP_e.N_mp]
            vzn1 = MP_e.vz_mp[0:MP_e.N_mp]

            if Ez_n != 0.:
                raise ValueError('Oooops! Not implemented....')

            boris_step_multipole(N_sub_steps, Dt_substep, self.B_field, self.B_field_skew,
                         xn1, yn1, zn1, vxn1, vyn1, vzn1,
                         Ex_n, Ey_n, Bx_arr, By_arr, Bz_arr, custom_B, MP_e.charge, MP_e.mass)

        return MP_e

#----------------------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 6.7.0
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

from numpy import sqrt, exp
from numpy.random import rand
import numpy as np

def yield_fun2(E, costheta, Emax, del_max, R0, E0, s, flag_costheta_delta_scale=True, flag_costheta_Emax_shift=True):

    if flag_costheta_delta_scale:
        del_max_tilde=del_max*exp(0.5*(1.-costheta))
    else:
        del_max_tilde=del_max

    if flag_costheta_Emax_shift:
        E_max_tilde=Emax*(1.+0.7*(1.-costheta))
    else:
        E_max_tilde=Emax

    x=E/E_max_tilde

    true_sec=del_max_tilde*(s*x)/(s-1.+x**s)
    reflected=R0*((sqrt(E)-sqrt(E+E0))/(sqrt(E)+sqrt(E+E0)))**2.

    delta = np.ones_like(E, float)
    flag_true_sec = reflected < rand(E.size)
    delta[flag_true_sec] = true_sec[flag_true_sec] / (1.-reflected[flag_true_sec])


    #print(delta.sum(), np.sum(true_sec+reflected))

    return delta, ~flag_true_sec


class SEY_model_ECLOUD:
    def __init__(self, Emax,del_max,R0,E0=150., s=1.35, flag_costheta_delta_scale=True, flag_costheta_Emax_shift=True):
            self.Emax = Emax
            self.del_max = del_max
            self.R0 = R0
            self.E0 = E0
            self.s = s
            self.flag_costheta_delta_scale = flag_costheta_delta_scale
            self.flag_costheta_Emax_shift = flag_costheta_Emax_shift

            print 'Secondary emission model: ECLOUD E0=%.4f s=%.4f' % (self.E0, self.s)

    def SEY_process(self, nel_impact, E_impact_eV, costheta_impact, i_impact):
            yiel, flag_elast=yield_fun2(E_impact_eV,costheta_impact,self.Emax,self.del_max,self.R0, E0=self.E0, s=self.s,
                                        flag_costheta_delta_scale=self.flag_costheta_delta_scale,
                                        flag_costheta_Emax_shift=self.flag_costheta_Emax_shift)
            nel_emit=nel_impact*yiel

            return nel_emit, flag_elast, ~flag_elast


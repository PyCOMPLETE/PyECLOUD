#----------------------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 6.4.0
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

from numpy import sqrt, exp, take
from numpy.random import rand
import numpy as np

def yield_fun2(E,costheta,Emax,del_max,R0,E0):

    s=1.35;

    del_max_tilde=del_max*exp(0.5*(1.-costheta));
    E_max_tilde=Emax*(1.+0.7*(1.-costheta));

    x=E/E_max_tilde;

    true_sec=del_max_tilde*(s*x)/(s-1.+x**s);
    reflected=R0*((sqrt(E)-sqrt(E+E0))/(sqrt(E)+sqrt(E+E0)))**2.;

    delta=true_sec+reflected;
    ref_frac=0.*delta
    mask_non_zero=(delta>0)
    ref_frac[mask_non_zero]=reflected[mask_non_zero]/delta[mask_non_zero];

    return delta, ref_frac


class SEY_model_ECLOUD_non_unif:
    def __init__(self, chamb, Emax,del_max,R0,E0=150.):
            if chamb.chamb_type!='polyg':
                raise ValueError("""ECLOUD_nunif can be used only with chamb_type='polyg'!!!""")

            self.del_max_segments = np.float_(chamb.del_max_segments)
            self.R0_segments = np.float_(chamb.R0_segments)
            self.Emax_segments = np.float_(chamb.Emax_segments)

            self.del_max_segments[chamb.del_max_segments<0.]=del_max
            self.R0_segments[chamb.R0_segments<0.]=R0
            self.Emax_segments[chamb.Emax_segments<0.]=Emax

            self.E0=E0

            print('Secondary emission model: ECLOUD non uniform E0=%f'%self.E0)

    def SEY_process(self,nel_impact,E_impact_eV, costheta_impact, i_impact):


            Emax_mp = take(self.Emax_segments, i_impact)
            del_max_mp = take(self.del_max_segments, i_impact)
            R0_mp = take(self.R0_segments, i_impact)


            yiel, ref_frac=yield_fun2(E_impact_eV,costheta_impact,Emax_mp,del_max_mp,R0_mp, E0=self.E0);
            flag_elast=(rand(len(ref_frac))<ref_frac);
            flag_truesec=~(flag_elast);
            nel_emit=nel_impact*yiel;

            return  nel_emit, flag_elast, flag_truesec


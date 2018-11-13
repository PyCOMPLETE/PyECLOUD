#-Begin-preamble-------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 7.6.0
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
#                           Lotta Mether
#                           Annalisa Romano
#                           Giovanni Rumolo
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

from numpy import sqrt, exp
from numpy.random import rand
from sec_emission_model_ECLOUD import SEY_model_ECLOUD


def yield_fun_furman_pivi(E,costheta,Emax,del_max,R0,E0):
    del_true_sec = None
    del_reflected = None
    del_rediffused = None
    return del_true_sec+del_reflected+del_rediffused, del_reflected, del_rediffused


class SEY_model_furman_pivi(SEY_model_ECLOUD):
    def __init__(
                    self, Emax, del_max, R0,
                    E_th=None, sigmafit=None, mufit=None,
                    switch_no_increase_energy=0, thresh_low_energy=None, secondary_angle_distribution=None,
                    E0=150., s=1.35, flag_costheta_delta_scale=True, flag_costheta_Emax_shift=True
                ):

        self.E_th = E_th
        self.sigmafit = sigmafit
        self.mufit = mufit
        self.switch_no_increase_energy = switch_no_increase_energy
        self.thresh_low_energy = thresh_low_energy
        self.secondary_angle_distribution = secondary_angle_distribution

        if secondary_angle_distribution is not None:
            import electron_emission
            self.angle_dist_func = electron_emission.get_angle_dist_func(secondary_angle_distribution)
        else:
            self.angle_dist_func = None

        self.Emax = Emax
        self.del_max = del_max
        self.R0 = R0
        self.E0 = E0
        self.s = s
        self.flag_costheta_delta_scale = flag_costheta_delta_scale
        self.flag_costheta_Emax_shift = flag_costheta_Emax_shift

        print 'Secondary emission model: Furman-Pivi E0=%.4f s=%.4f' % (self.E0, self.s)

    def SEY_process(self,nel_impact,E_impact_eV, costheta_impact, i_impact):

        yiel, ref_frac, del_rediffused = yield_fun_furman_pivi(
            E_impact_eV, costheta_impact,self.Emax, self.del_max, self.R0, E0=self.E0, s=self.s,
            flag_costheta_delta_scale=self.flag_costheta_delta_scale, flag_costheta_Emax_shift=self.flag_costheta_Emax_shift)
        flag_elast = None
        flag_truesec = None
        nel_emit = None

        return nel_emit, flag_elast, flag_truesec

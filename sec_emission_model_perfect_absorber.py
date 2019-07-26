#-Begin-preamble-------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 8.1.0
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


from __future__ import division, print_function
import os

import numpy as np
import scipy.io as sio
from numpy.random import rand
from sec_emission_model_ECLOUD import SEY_model_ECLOUD


class SEY_model_perfect_absorber(object):
    event_types = {
        0: 'absorb',
    }


    def __init__(self):
        """
        Perfect absorber
        """
        pass

    def SEY_model_evol(self, Dt):
        pass

    def impacts_on_surface(self, mass, nel_impact, x_impact, y_impact, z_impact,
                           vx_impact, vy_impact, vz_impact, Norm_x, Norm_y, i_found,
                           v_impact_n, E_impact_eV, costheta_impact, nel_mp_th, flag_seg):

        nel_replace = nel_impact * 0.
        x_replace = x_impact.copy()
        y_replace = y_impact.copy()
        z_replace = z_impact.copy()
        vx_replace = vx_impact.copy()
        vy_replace = vy_impact.copy()
        vz_replace = vz_impact.copy()

        if i_found is not None:
            i_seg_replace = i_found.copy()
        else:
            i_seg_replace = i_found

        
        nel_new_MPs = np.array([])
        x_new_MPs = np.array([])
        y_new_MPs = np.array([])
        z_new_MPs = np.array([])
        vx_new_MPs = np.array([])
        vy_new_MPs = np.array([])
        vz_new_MPs = np.array([])
        i_seg_new_MPs = np.array([])

        nel_emit_tot_events = nel_impact * 0.
        event_type = np.zeros_like(nel_impact)

        # extended_event_type keeps track of the event type for new MPs, it is
        # needed for the extraction of emission-energy distributions.

        event_info = {'extended_event_type': event_type}

        return nel_emit_tot_events, event_type, event_info,\
            nel_replace, x_replace, y_replace, z_replace, vx_replace, vy_replace, vz_replace, i_seg_replace,\
            nel_new_MPs, x_new_MPs, y_new_MPs, z_new_MPs, vx_new_MPs, vy_new_MPs, vz_new_MPs, i_seg_new_MPs




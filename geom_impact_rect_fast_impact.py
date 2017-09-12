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


import geom_impact_poly_fast_impact as gipfi
import numpy as np

na = np.array


def rect_cham_geom_object(x_aper, y_aper, flag_verbose_file=False, flag_verbose_stdout=False):
    chamber = gipfi.polyg_cham_geom_object(
        {'Vx':na([x_aper, -x_aper, -x_aper, x_aper]),
            'Vy':na([y_aper, y_aper, -y_aper, -y_aper]),
            'x_sem_ellip_insc':0.99*x_aper,
            'y_sem_ellip_insc':0.99*y_aper},
        flag_non_unif_sey = False,
        flag_verbose_file=flag_verbose_file, flag_verbose_stdout=flag_verbose_stdout)

    chamber.chamb_type='rect'

    return chamber


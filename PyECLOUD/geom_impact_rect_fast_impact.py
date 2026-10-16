# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

from . import geom_impact_poly_fast_impact as gipfi
import numpy as np

na = np.array


def rect_cham_geom_object(x_aper, y_aper, flag_non_unif_sey, **kwargs):
    chamber = gipfi.polyg_cham_geom_object(
        {
            'Vx': na([x_aper, -x_aper, -x_aper, x_aper]),
            'Vy': na([y_aper, y_aper, -y_aper, -y_aper]),
            'x_sem_ellip_insc': 0.99 * x_aper,
            'y_sem_ellip_insc': 0.99 * y_aper
        }, flag_non_unif_sey, **kwargs)

    chamber.chamb_type = 'rect'

    return chamber


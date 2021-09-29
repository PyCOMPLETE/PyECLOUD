

import sys
sys.path.append("../")

import numpy as np
from scipy.constants import c, e, m_p, m_e

from PyHEADTAIL.machines.synchrotron import Synchrotron


class EmptyObject(object):
    pass


class CLIC_DR(Synchrotron):

    def __init__(self, machine_configuration=None, optics_mode='smooth', wrap_z=True, **kwargs):


        pp = EmptyObject()

        pp.machine_configuration = machine_configuration
        pp.optics_mode = optics_mode
        pp.wrap_z = wrap_z

        pp.longitudinal_mode = 'non-linear' # we have to take the non-linear to try the multibunch
        pp.alpha        = 0.0001276102729
        pp.p0           = 2.86e9 * e / c
        pp.p_increment = 0.
        pp.accQ_x       = 48.35
        pp.accQ_y       = 10.40
        pp.mass         = m_e
        pp.charge       = -e
        pp.RF_at        = 'end_of_transverse'

        if machine_configuration == 'CLIC_DR_1GHz':
            pp.h_RF     = 1425
            pp.V_RF = 5.1e6
            pp.dphi_RF = 0

        elif machine_configuration == 'CLIC_DR_2GHz':
            pp.h_RF     = 2851
            pp.V_RF = 4.5e6
            pp.dphi_RF = 0

        else:
            raise ValueError('ERROR: unknown machine configuration', pp.machine_configuration)

        if optics_mode == 'smooth':
            if 's' in list(kwargs.keys()):
                raise ValueError('s vector cannot be provided if optics_mode = "smooth"')

            pp.n_segments = kwargs['n_segments']
            pp.circumference = 427.5

            pp.name = None

            pp.beta_x       = 4.0
            pp.D_x      = 0
            pp.beta_y       = 9.0
            pp.D_y      = 0

            pp.alpha_x  = None
            pp.alpha_y  = None

            pp.s = None

        elif optics_mode == 'non-smooth':
            if 'n_segments' in list(kwargs.keys()):
                raise ValueError('n_segments cannot be provided if optics_mode = "non-smooth"')
            pp.n_segments = None
            pp.circumference = None

            try:
                pp.name     = kwargs['name']
            except KeyError:
                pp.name = None

            pp.beta_x       = kwargs['beta_x']
            pp.beta_y       = kwargs['beta_y']

            try:
                pp.D_x  = kwargs['D_x']
            except KeyError:
                pp.D_x  = 0 * np.array(kwargs['s'])
            try:
                pp.D_y  = kwargs['D_y']
            except KeyError:
                pp.D_y  = 0 * np.array(kwargs['s'])

            pp.alpha_x  = kwargs['alpha_x']
            pp.alpha_y  = kwargs['alpha_y']

            pp.accQ_x       = kwargs['accQ_x']
            pp.accQ_y       = kwargs['accQ_y']

            pp.s = kwargs['s']

        else:
            raise ValueError('optics_mode not recognized!')

        # detunings
        pp.Qp_x     = 0.
        pp.Qp_y     = 0.

        pp.app_x        = 0
        pp.app_y        = 0
        pp.app_xy       = 0

        for attr in list(kwargs.keys()):
            if kwargs[attr] is not None:
                if type(kwargs[attr]) is list or type(kwargs[attr]) is np.ndarray:
                    str2print = '[%s ...]'%repr(kwargs[attr][0])
                else:
                    str2print = repr(kwargs[attr])
                self.prints('Synchrotron init. From kwargs: %s = %s'
                                    % (attr, str2print))

                if not hasattr(pp, attr):
                    raise NameError("I don't understand %s"%attr)

                setattr(pp, attr, kwargs[attr]) 

        super(CLIC_DR, self).__init__(optics_mode=pp.optics_mode, circumference=pp.circumference, n_segments=pp.n_segments, s=pp.s, name=pp.name,
                                alpha_x=pp.alpha_x, beta_x=pp.beta_x, D_x=pp.D_x, alpha_y=pp.alpha_y, beta_y=pp.beta_y, D_y=pp.D_y,
                                accQ_x=pp.accQ_x, accQ_y=pp.accQ_y, Qp_x=pp.Qp_x, Qp_y=pp.Qp_y, app_x=pp.app_x, app_y=pp.app_y, app_xy=pp.app_xy,
                                alpha_mom_compaction=pp.alpha, longitudinal_mode=pp.longitudinal_mode,
                                h_RF=np.atleast_1d(pp.h_RF), V_RF=np.atleast_1d(pp.V_RF), dphi_RF=np.atleast_1d(pp.dphi_RF), p0=pp.p0, p_increment=pp.p_increment,
                                charge=pp.charge, mass=pp.mass, RF_at=pp.RF_at, wrap_z=pp.wrap_z)

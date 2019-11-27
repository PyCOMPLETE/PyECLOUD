from PyHEADTAIL.machines.synchrotron import Synchrotron
import numpy as np
from scipy.constants import c, e, m_p

class EmptyObject(object):
    pass

class LHC(Synchrotron):

    def __init__(self, machine_configuration=None, optics_mode='smooth', **kwargs):

        pp = EmptyObject()

        pp.machine_configuration = machine_configuration
        pp.optics_mode = optics_mode
        
        pp.longitudinal_mode = 'non-linear'
        pp.alpha       = 3.225e-04
        pp.h_RF        = 35640
        pp.mass        = m_p
        pp.charge      = e

        if pp.machine_configuration == 'Injection':
            pp.p0          = 450e9 * e / c
            pp.p_increment = 0.
            pp.accQ_x      = 64.28
            pp.accQ_y      = 59.31
            pp.V_RF        = 6e6
            pp.dphi_RF     = 0.
        elif pp.machine_configuration == '6.5_TeV_collision_tunes':
            pp.p0          = 6500e9 * e / c
            pp.p_increment = 0.
            pp.accQ_x      = 64.31
            pp.accQ_y      = 59.32
            pp.V_RF        = 12e6
            pp.dphi_RF     = 0.
        else:
            raise ValueError('machine_configuration not recognized!')

        if pp.optics_mode == 'smooth':
            if 's' in list(kwargs.keys()):
                raise ValueError('s vector cannot be provided if optics_mode = "smooth"')

            pp.n_segments = kwargs['n_segments']
            pp.circumference = 26658.8832

            pp.name = None

            pp.beta_x      = 92.7
            pp.D_x         = 0
            pp.beta_y      = 93.2
            pp.D_y         = 0

            pp.alpha_x     = None
            pp.alpha_y     = None

            pp.s = None

        elif pp.optics_mode == 'non-smooth':
            if 'n_segments' in list(kwargs.keys()):
                raise ValueError('n_segments cannot be provided if optics_mode = "non-smooth"')
            pp.n_segments = None
            pp.circumference = None

            pp.name        = kwargs['name']

            pp.beta_x      = kwargs['beta_x']
            pp.beta_y      = kwargs['beta_y']

            try:
                pp.D_x         = kwargs['D_x']
            except KeyError:
                pp.D_x         = 0 * np.array(kwargs['s'])
            try:
                pp.D_y         = kwargs['D_y']
            except KeyError:
                pp.D_y         = 0 * np.array(kwargs['s'])

            pp.alpha_x     = kwargs['alpha_x']
            pp.alpha_y     = kwargs['alpha_y']

            pp.s = kwargs['s']

        else:
            raise ValueError('optics_mode not recognized!')

        # detunings
        pp.Qp_x        = 0
        pp.Qp_y        = 0

        pp.app_x       = 0
        pp.app_y       = 0
        pp.app_xy      = 0

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


        super(LHC, self).__init__(optics_mode=pp.optics_mode, circumference=pp.circumference, n_segments=pp.n_segments, s=pp.s, name=pp.name,
                            alpha_x=pp.alpha_x, beta_x=pp.beta_x, D_x=pp.D_x, alpha_y=pp.alpha_y, beta_y=pp.beta_y, D_y=pp.D_y,
                            accQ_x=pp.accQ_x, accQ_y=pp.accQ_y, Qp_x=pp.Qp_x, Qp_y=pp.Qp_y, app_x=pp.app_x, app_y=pp.app_y, app_xy=pp.app_xy,
                            alpha_mom_compaction=pp.alpha, longitudinal_mode=pp.longitudinal_mode,
                            h_RF=np.atleast_1d(pp.h_RF), V_RF=np.atleast_1d(pp.V_RF), dphi_RF=np.atleast_1d(pp.dphi_RF), p0=pp.p0, p_increment=pp.p_increment,
                            charge=pp.charge, mass=pp.mass)



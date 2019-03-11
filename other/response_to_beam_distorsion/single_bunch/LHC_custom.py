from PyHEADTAIL.machines.synchrotron import BasicSynchrotron
import numpy as np
from scipy.constants import c, e, m_p


class EmptyObject(object):
    pass


class LHC(BasicSynchrotron):

    def __init__(self, n_segments, machine_configuration, **kwargs):


        pp = EmptyObject()

        pp.machine_configuration = machine_configuration

        pp.circumference = 35640 * 2.5e-9 * c
        pp.longitudinal_mode = 'non-linear'
        pp.p_increment = 0.
        pp.charge = e
        pp.mass = m_p
        pp.alpha = 3.225e-04

        if machine_configuration == 'HLLHC-injection':
            pp.alpha_x = 0.
            pp.beta_x = 92.7
            pp.D_x = 0.
            pp.alpha_y = 0.
            pp.beta_y = 93.2
            pp.D_y = 0.

            pp.accQ_x = 62.28
            pp.accQ_y = 60.31

            pp.h_RF = 35640
            pp.V_RF = 8e6
            pp.dphi_RF = 0.

            pp.p0 = 450.e9 * e / c

        elif machine_configuration == 'HLLHC-collision':
            pp.alpha_x = 0.
            pp.beta_x = 92.7
            pp.D_x = 0.
            pp.alpha_y = 0.
            pp.beta_y = 93.2
            pp.D_y = 0.

            pp.accQ_x = 62.31
            pp.accQ_y = 60.32

            pp.h_RF = 35640
            pp.V_RF = 16e6
            pp.dphi_RF = 0.

            pp.p0 = 7000e9 * e / c

        else:
            raise ValueError('ERROR: unknown machine configuration', pp.machine_configuration)

        # detunings
        pp.Qp_x = 0
        pp.Qp_y = 0

        pp.app_x = 0
        pp.app_y = 0
        pp.app_xy = 0

        for attr in kwargs.keys():
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


        super(LHC, self).__init__(optics_mode='smooth', circumference=pp.circumference, n_segments=pp.n_segments,
                                  alpha_x=pp.alpha_x, beta_x=pp.beta_x, D_x=pp.D_x, alpha_y=pp.alpha_y, beta_y=pp.beta_y, D_y=pp.D_y,
                                  accQ_x=pp.accQ_x, accQ_y=pp.accQ_y, Qp_x=pp.Qp_x, Qp_y=pp.Qp_y, app_x=pp.app_x, app_y=pp.app_y, app_xy=pp.app_xy,
                                  alpha_mom_compaction=pp.alpha, longitudinal_mode=pp.longitudinal_mode,
                                  h_RF=np.atleast_1d(pp.h_RF), V_RF=np.atleast_1d(pp.V_RF), dphi_RF=np.atleast_1d(pp.dphi_RF), p0=pp.p0, p_increment=pp.p_increment,
                                  charge=pp.charge, mass=pp.mass, RF_at='end_of_transverse')


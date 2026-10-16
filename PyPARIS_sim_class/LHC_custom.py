# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

from PyHEADTAIL.machines.synchrotron import BasicSynchrotron
import numpy as np
from scipy.constants import c, e, m_p

class EmptyObject(object):
    pass

class LHC(BasicSynchrotron):

    def __init__(self, n_segments, machine_configuration, optics_dict=None, **kwargs):
        
        
        pp = EmptyObject()

        pp.n_segments = n_segments
        pp.machine_configuration = machine_configuration
        pp.optics_dict = optics_dict

        pp.circumference     = 26658.8832
        pp.longitudinal_mode = 'non-linear'
        pp.p_increment       = 0.
        pp.charge            = e
        pp.mass              = m_p
        pp.alpha             = 3.225e-04
        

        if pp.machine_configuration =='HLLHC-injection':
            pp.alpha_x        = 0. 
            pp.beta_x         = 92.7 
            pp.D_x            = 0. 
            pp.alpha_y        = 0. 
            pp.beta_y         = 93.2 
            pp.D_y            = 0.

            pp.accQ_x         = 62.27
            pp.accQ_y         = 60.295
            
            pp.h_RF           = 35640
            pp.V_RF           = 8e6
            pp.dphi_RF        = 0.
            
            pp.p0 = 450.e9 * e /c

        elif pp.machine_configuration =='HLLHC-collision':
            pp.alpha_x        = 0. 
            pp.beta_x         = 92.7 
            pp.D_x            = 0. 
            pp.alpha_y        = 0. 
            pp.beta_y         = 93.2 
            pp.D_y            = 0.

            pp.accQ_x         = 62.31
            pp.accQ_y         = 60.32
            
            pp.h_RF           = 35640
            pp.V_RF           = 16e6
            pp.dphi_RF        = 0.
            
            pp.p0 = 7000e9 * e /c
            
        elif machine_configuration =='LHC-collision':
            pp.alpha_x        = 0. 
            pp.beta_x         = 92.7 
            pp.D_x            = 0. 
            pp.alpha_y        = 0. 
            pp.beta_y         = 93.2 
            pp.D_y            = 0.

            pp.accQ_x         = 62.31
            pp.accQ_y         = 60.32
            
            pp.h_RF           = 35640
            pp.V_RF           = 12e6
            pp.dphi_RF        = 0.
            
            pp.p0 = 6800e9 * e /c
            

        else:
            raise ValueError('ERROR: unknown machine configuration', machine_configuration)

        # detunings
        pp.Qp_x        = 0
        pp.Qp_y        = 0
        
        pp.app_x       = 0
        pp.app_y       = 0
        pp.app_xy      = 0
        
        pp.i_octupole_focusing = None
        pp.i_octupole_defocusing = None
        pp.octupole_knob = None	
        
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
                

        if pp.i_octupole_focusing is not None or pp.i_octupole_defocusing is not None:
            if pp.octupole_knob is not None:
                raise ValueError('octupole_knobs and octupole currents cannot be used at the same time!')
            pp.app_x, pp.app_y, pp.app_xy = self._anharmonicities_from_octupole_current_settings(
                pp.i_octupole_focusing, pp.i_octupole_defocusing)
            self.i_octupole_focusing = pp.i_octupole_focusing
            self.i_octupole_defocusing = pp.i_octupole_defocusing
            
        if pp.octupole_knob is not None:	
            if pp.i_octupole_focusing is not None or pp.i_octupole_defocusing is not None:
                raise ValueError('octupole_knobs and octupole currents cannot be used at the same time!')
            pp.i_octupole_focusing, pp.i_octupole_defocusing = self._octupole_currents_from_octupole_knobs(pp.octupole_knob, pp.p0)
            pp.app_x, pp.app_y, pp.app_xy = self._anharmonicities_from_octupole_current_settings(
                pp.i_octupole_focusing, pp.i_octupole_defocusing)		
            self.i_octupole_focusing = pp.i_octupole_focusing
            self.i_octupole_defocusing = pp.i_octupole_defocusing
            
        if pp.optics_dict is not None:
            if pp.n_segments is not None: raise ValueError('n_segments cannot be provided if optics_mode = "non-smooth"')
            
            for vv in 'beta_x beta_y D_x D_y alpha_x alpha_y s accQ_x accQ_y'.split():
                if vv in kwargs:
                    if kwargs[vv] is not None:
                        raise ValueError('%s cannot be provided when optics_dict is given!'%vv)
            

            pp.n_segments = None
            pp.circumference = None
            pp.name = pp.optics_dict['name']
            
            pp.beta_x = pp.optics_dict['beta_x']
            pp.beta_y = pp.optics_dict['beta_y'] 
            
            try:
                pp.D_x = pp.optics_dict['D_x']	
            except KeyError:
                pp.D_x = 0*np.array(pp.optics_dict['s'])
            try:
                pp.D_y = kwargs['D_y']	
            except KeyError:
                pp.D_y = 0*np.array(pp.optics_dict['s'])
                    
            pp.alpha_x = pp.optics_dict['alpha_x']
            pp.alpha_y = pp.optics_dict['alpha_y']
            
            pp.s = pp.optics_dict['s']
            
            pp.accQ_x = pp.optics_dict['accQ_x']
            pp.accQ_y = pp.optics_dict['accQ_y']
            pp.optics_mode = 'non-smooth'
        else:
            pp.optics_mode = 'smooth'
            pp.s = None
            pp.name = None
        
        
        
        
        super(LHC, self).__init__(optics_mode=pp.optics_mode, circumference=pp.circumference, n_segments=pp.n_segments,
             s=pp.s, name=pp.name, 
             alpha_x=pp.alpha_x, beta_x=pp.beta_x, D_x=pp.D_x, alpha_y=pp.alpha_y, beta_y=pp.beta_y, D_y=pp.D_y,
             accQ_x=pp.accQ_x, accQ_y=pp.accQ_y, Qp_x=pp.Qp_x, Qp_y=pp.Qp_y, app_x=pp.app_x, app_y=pp.app_y, app_xy=pp.app_xy,
             alpha_mom_compaction=pp.alpha, longitudinal_mode=pp.longitudinal_mode,
             h_RF=np.atleast_1d(pp.h_RF), V_RF=np.atleast_1d(pp.V_RF), dphi_RF=np.atleast_1d(pp.dphi_RF), 
             p0=pp.p0, p_increment=pp.p_increment,
             charge=pp.charge, mass=pp.mass, RF_at='end_of_transverse')
             
    def _anharmonicities_from_octupole_current_settings(self, i_octupole_focusing, i_octupole_defocusing):
                """Calculate the constants of proportionality app_x, app_y and
                app_xy (== app_yx) for the amplitude detuning introduced by the
                LHC octupole magnets (aka. LHC Landau octupoles) from the
                electric currents i_octupole_focusing [A] and i_octupole_defocusing [A] flowing
                through the magnets. The maximum current is given by
                i_max = +/- 550 [A]. The values app_x, app_y, app_xy obtained
                from the formulae are proportional to the strength of detuning
                for one complete turn around the accelerator, i.e. one-turn
                values.
                The calculation is based on formulae (3.6) taken from 'The LHC
                transverse coupled-bunch instability' by N. Mounet, EPFL PhD
                Thesis, 2012. Values (hard-coded numbers below) are valid for
                LHC Landau octupoles before LS1. Beta functions in x and y are
                correctly taken into account. Note that here, the values of
                app_x, app_y and app_xy are not normalized to the reference
                momentum p0. This is done only during the calculation of the
                detuning in the corresponding detune method of the
                AmplitudeDetuningSegment.
                More detailed explanations and references on how the formulae
                were obtained are given in the PhD thesis (pg. 85ff) cited
                above.
                """
                i_max = 550.  # [A]
                E_max = 7000. # [GeV]
        
                app_x  = E_max * (267065. * i_octupole_focusing / i_max -
                    7856. * i_octupole_defocusing / i_max)
                app_y  = E_max * (9789. * i_octupole_focusing / i_max -
                    277203. * i_octupole_defocusing / i_max)
                app_xy = E_max * (-102261. * i_octupole_focusing / i_max +
                    93331. * i_octupole_defocusing / i_max)
        
                # Convert to SI units.
                convert_to_SI = e / (1.e-9 * c)
                app_x *= convert_to_SI
                app_y *= convert_to_SI
                app_xy *= convert_to_SI
        
                return app_x, app_y, app_xy
            
    def _octupole_currents_from_octupole_knobs(self, octupole_knob, p0):
        i_octupole_focusing = 19.557 * octupole_knob / (-1.5) * p0 / 2.4049285931335872e-16
        i_octupole_defocusing = - i_octupole_focusing
        return i_octupole_focusing, i_octupole_defocusing


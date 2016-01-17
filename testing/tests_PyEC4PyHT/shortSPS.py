from PyCERNmachines.machines import synchrotron
import numpy as np
from scipy.constants import c, e, m_p

class shortSPS(synchrotron):

	def __init__(self, *args, **kwargs):
		
		if 'n_segments' not in kwargs.keys():
			raise ValueError('Number of segments must be specified')
			
		if 'machine_configuration' not in kwargs.keys():
			raise ValueError('machine_configuration must be specified')
			
		self.n_segments = kwargs['n_segments']
		self.machine_configuration = kwargs['machine_configuration']
		
		self.circumference = 6911./6.

		self.s = np.arange(0, self.n_segments + 1) * self.circumference / self.n_segments

		if self.machine_configuration =='Q20-injection-like':
			self.charge = e
			self.mass = m_p
			
			self.alpha_x        = 0 * np.ones(self.n_segments)
			self.beta_x         = 54.6 * np.ones(self.n_segments)
			self.D_x            = 0 * np.ones(self.n_segments)
			self.alpha_y        = 0 * np.ones(self.n_segments)
			self.beta_y         = 54.6 * np.ones(self.n_segments)
			self.D_y            = 0 * np.ones(self.n_segments)

			self.Q_x            = 5.13
			self.Q_y            = 5.18

			self.Qp_x           = 0
			self.Qp_y           = 0

			self.app_x          = 0.0000e-9
			self.app_y          = 0.0000e-9
			self.app_xy         = 0

			self.alpha     		= 6*0.00308642
			
			self.h1, self.h2       = 4620/6, 4620*4/6
			self.V1, self.V2       = 5.75e6, 0
			self.dphi1, self.dphi2 = 0, np.pi
			
			
			self.longitudinal_focusing = 'linear'
			
			self.gamma = 27.7
			
			self.p_increment       = 0 * e/c * self.circumference/(self.beta*c)
			
		else:
			raise ValueError('ERROR: unknown machine configuration', machine_configuration)

		
		super(shortSPS, self).__init__(*args, **kwargs)

# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

import numpy as np

class ShareSegments(object):
	def __init__(self, N_segments, N_nodes):
		self.N_segments_per_node = np.array(N_nodes*[N_segments//N_nodes])
		rem = N_segments-N_segments//N_nodes*N_nodes
		self.N_segments_per_node[:rem]+=1
		self.start_parts = np.array([0]+list(np.cumsum(self.N_segments_per_node)))
		
	def my_part(self, myid):
		return self.start_parts[myid], self.start_parts[myid+1]


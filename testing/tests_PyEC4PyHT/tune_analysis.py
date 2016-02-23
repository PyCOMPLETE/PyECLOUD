import numpy as np
from PySUSSIX import Sussix
import sys

def tune_analysis(x_i, xp_i, y_i, yp_i):
		n_turns = x_i.shape[1]
		macroparticlenumber = x_i.shape[0]
		
		qx_i = np.empty(macroparticlenumber)
		qy_i = np.empty(macroparticlenumber)

		spectrum_x = np.abs(np.fft.fft(x_i, axis = 1))[:,:n_turns/2]
		tune_peak_x = np.float_(np.argmax(spectrum_x, axis=1))/float(n_turns)

		spectrum_y = np.abs(np.fft.fft(y_i, axis = 1))[:,:n_turns/2]
		tune_peak_y = np.float_(np.argmax(spectrum_y, axis=1))/float(n_turns)

		print 'analysing particle spectra ... this may take some time.'
		for p_idx in xrange(macroparticlenumber):
			
			SX = Sussix()
			SX.sussix_inp(nt1=1, nt2=n_turns, idam=2, ir=0, tunex=tune_peak_x[p_idx], tuney=tune_peak_y[p_idx])

			SX.sussix(x_i[p_idx,:], xp_i[p_idx,:],
					  y_i[p_idx,:], yp_i[p_idx,:],
					  x_i[p_idx,:], xp_i[p_idx,:]) 
			qx_i[p_idx] = SX.ox[0]
			qy_i[p_idx] = SX.oy[0]
			
			sys.stdout.write('\rparticle %d'%p_idx)
				
		x_centroid = np.mean(x_i, axis=0)
		xp_centroid = np.mean(xp_i, axis=0)
		y_centroid = np.mean(y_i, axis=0)
		yp_centroid = np.mean(yp_i, axis=0)
		
		#print x_centroid.shape
		
		spectrum_x_centroid = np.abs(np.fft.fft(x_centroid))[:n_turns/2]
		tune_peak_x_centroid = np.float_(np.argmax(spectrum_x_centroid))/float(n_turns)

		spectrumy_centroid = np.abs(np.fft.fft(x_centroid))[:n_turns/2]
		tune_peak_y_centroid = np.float_(np.argmax(spectrumy_centroid))/float(n_turns)

		SX = Sussix()
		SX.sussix_inp(nt1=1, nt2=n_turns, idam=2, ir=0, tunex=tune_peak_x_centroid, tuney=tune_peak_y_centroid)
		SX.sussix(x_centroid, xp_centroid,
				  y_centroid, yp_centroid,
				  x_centroid, xp_centroid) 
		qx_centroid = SX.ox[0]
		qy_centroid = SX.oy[0]

		return qx_i, qy_i, qx_centroid, qy_centroid 



# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

import numpy as np
from PyHEADTAIL.particles.particles import Particles
import json
import pickle


def combine_float_buffers(list_of_buffers):
	N_buffers = len(list_of_buffers)
	len_buffers = np.array([float(len(seq)) for seq in list_of_buffers])
	return np.array(np.concatenate([np.array([float(N_buffers)]), len_buffers]+list_of_buffers),dtype=np.float64)

def split_float_buffers(megabuffer):
	i_mbuf = 0
	
	N_buffers = int(megabuffer[0])
	i_mbuf += 1
	
	len_buffers = megabuffer[i_mbuf:i_mbuf+N_buffers]
	i_mbuf += N_buffers
	
	list_of_buffers = []
	for i_buf in range(N_buffers):
		lenbuf = int(len_buffers[i_buf])
		list_of_buffers.append(megabuffer[i_mbuf:i_mbuf+lenbuf])
		i_mbuf += lenbuf
		
	return list_of_buffers
	
	
	


def list_of_strings_2_buffer(strlist):
	data = ''.join([s+';' for s in strlist])+'\nendbuf\n'
	buf_to_send = np.atleast_1d(np.int_(np.array(list(map(ord, list(data))))))
	return buf_to_send
	
def buffer_2_list_of_strings(buf):
	str_received = ''.join(map(chr, list(buf)))
	strlist = list(map(str, str_received.split('\nendbuf\n')[0].split(';')))[:-1]
	return strlist


def beam_2_buffer(beam, mode='pickle', verbose=False):
	
	#print beam
	
	if beam is None:
		buf = np.array([-1.])
	else:	
		if np.array(beam.particlenumber_per_mp).shape != ():
			raise ValueError('particlenumber_per_mp is a vector! Not implemented!')

		
		# if not hasattr(beam, 'slice_info'):
		# 	sl_info_buf = np.array([-1., 0., 0., 0.])
		# elif beam.slice_info is None:
		# 	sl_info_buf = np.array([-1., 0., 0., 0.])
		# elif beam.slice_info == 'unsliced':
		# 	sl_info_buf = np.array([0., 0., 0., 0.])
		# else:
		# 	sl_info_buf = np.array([1.,
		# 					beam.slice_info['z_bin_center'],
		# 					beam.slice_info['z_bin_right'],
		# 					beam.slice_info['z_bin_left']])

		# Prepare slice info buffer
		if not hasattr(beam, 'slice_info'):
			sinfo = None
		else:
			sinfo = beam.slice_info

		# Beam data buffer
		if mode=='json':
			sinfo_str = json.dumps(sinfo)
			sinfo_int = np.array(list(map(ord, sinfo_str)), dtype=np.int)
			sinfo_float_buf = sinfo_int.astype(np.float, casting='safe')
		elif mode=='pickle':
			pss = pickle.dumps(sinfo, protocol=2)
			# Pad to have a multiple of 8 bytes
			s1arr = np.frombuffer(pss, dtype='S1')
			ll = len(s1arr)
			s1arr_padded = np.concatenate((s1arr, np.zeros(8-ll%8, dtype='S1')))
			# Cast to array of floats
			f8arr = np.frombuffer(s1arr_padded, dtype=np.float64)
			sinfo_float_buf = np.concatenate((np.array([ll], dtype=np.float64),f8arr))
		else:
			raise ValueError('Unknown mode!')

		
		buf = np.array(np.concatenate((
			np.array([np.float64(beam.macroparticlenumber)]),
			np.array([np.float64(beam.particlenumber_per_mp)]), 
			np.array([beam.charge]),
			np.array([beam.mass]),
			np.array([beam.circumference]),
			np.array([beam.gamma]),
			np.atleast_1d(np.float64(beam.id)),
			beam.x, beam.xp, beam.y, beam.yp, beam.z, beam.dp,
			np.array([float(len(sinfo_float_buf))]),sinfo_float_buf)), dtype=np.float64)

		if verbose:
			print(('mode=%s'%mode))
			print(('beam.macroparticlenumber:%d'%beam.macroparticlenumber))
			print(('len(buf):%d'%len(buf)))
			print(('len(sinfo_float_buf):%d'%len(sinfo_float_buf)))
			
	return buf
	
def buffer_2_beam(buf, mode='pickle'):
	
	if buf[0]<0:
		beam=None
	else:
		i_buf = 0
		
		macroparticlenumber = int(buf[i_buf])
		i_buf += 1
		
		particlenumber_per_mp = buf[i_buf]
		i_buf += 1
		
		charge = buf[i_buf]
		i_buf += 1
		
		mass = buf[i_buf]
		i_buf += 1
		
		circumference = buf[i_buf]
		i_buf += 1
		
		gamma = buf[i_buf]
		i_buf += 1
		
		id_ = buf[i_buf:i_buf+macroparticlenumber].astype(np.int32)
		i_buf += macroparticlenumber

		x =  buf[i_buf:i_buf+macroparticlenumber]
		i_buf += macroparticlenumber
		
		xp =  buf[i_buf:i_buf+macroparticlenumber]
		i_buf += macroparticlenumber
		
		y =  buf[i_buf:i_buf+macroparticlenumber]
		i_buf += macroparticlenumber
		
		yp =  buf[i_buf:i_buf+macroparticlenumber]
		i_buf += macroparticlenumber
		
		z =  buf[i_buf:i_buf+macroparticlenumber]
		i_buf += macroparticlenumber	
		
		dp = buf[i_buf:i_buf+macroparticlenumber]
		i_buf += macroparticlenumber
		
		slice_info_buf_size = int(buf[i_buf])
		i_buf += 1
		
		slice_info_buf = buf[i_buf:i_buf+slice_info_buf_size]
		i_buf += slice_info_buf_size
		
		beam = Particles(macroparticlenumber=macroparticlenumber,
						particlenumber_per_mp=particlenumber_per_mp, charge=charge,
						mass=mass, circumference=circumference, gamma=gamma, 
						coords_n_momenta_dict={\
								'x': np.atleast_1d(x),
								'xp':np.atleast_1d(xp),
								'y':np.atleast_1d(y),
								'yp':np.atleast_1d(yp),	
								'z':np.atleast_1d(z),
								'dp':np.atleast_1d(dp)})
		
		beam.id = np.atleast_1d(id_)

		
		if mode=='json':
			si_int = slice_info_buf.astype(np.int)
			si_str = ''.join(map(chr, list(si_int)))
			beam.slice_info = json.loads(si_str)
		elif mode=='pickle':
			# Get length in bytes
			llrec = int(slice_info_buf[0])
			s1back_padded = np.frombuffer(slice_info_buf[1:].tobytes(), dtype='S1') 
			s1back = s1back_padded[:llrec]
			pss_rec = s1back.tobytes()
			beam.slice_info = pickle.loads(pss_rec)
		else:
			raise ValueError('Unknown mode!')


		
		# if slice_info_buf[0] < 0.:
		# 	beam.slice_info = None
		# elif slice_info_buf[0] == 0.:
		# 	beam.slice_info = 'unsliced'
		# else:
		# 	beam.slice_info = {\
  #                   'z_bin_center':slice_info_buf[1] ,
  #                   'z_bin_right':slice_info_buf[2],
  #                   'z_bin_left':slice_info_buf[3]}




	return beam

# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #


import scipy.io as sio
from pylab import imshow, colorbar, show, savefig, clf, squeeze, title, xlabel, ylabel, floor
import subprocess
import numpy as np

firs_pass = 0
last_pass = 3

flag_log = True

N_dec = 2

i_photog = 0

for pass_ind in range(firs_pass, last_pass):

    filename_rho = 'rho_video/rho_pass%d.mat'%pass_ind

    dict_ecl_video = sio.loadmat(filename_rho)
    dict_pyecltest = sio.loadmat('Pyecltest.mat')

    rho_video = dict_ecl_video['rho_video']
    t_video = squeeze(dict_ecl_video['t_video'].real)
    b_spac = squeeze(dict_pyecltest['b_spac'].real)
    (nphotog, _, _) = rho_video.shape

    #subprocess.check_call(('rm',  '*.png'))

    for ii in range(0, nphotog, N_dec):
        print('Pass %d %d/%d'%(pass_ind, ii, nphotog))
        imm = np.squeeze(rho_video[ii, :, :])
        if flag_log:
            imm = np.log10(np.abs(imm))
        imshow(imm.T, cmap=None, norm=None, aspect=None, interpolation=None,
               alpha=None, vmin=None, vmax=None, origin='lower', extent=None)
        colorbar()
        title(('passage = %d' % floor(t_video[ii] / b_spac)))
        #xlabel('x [m]')
        #ylabel('y [m]')
        filename = str('Pass%05d_%05d' % (pass_ind, ii)) + '.png'
        savefig(filename, dpi=100)
        clf()
        i_photog += 1


command = ('mencoder',
           'mf://*.png',
           '-mf',
           'type=png:w=800:h=600:fps=5',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4',
           '-oac',
           'copy',
           '-o',
           'output.avi')

subprocess.check_call(command)



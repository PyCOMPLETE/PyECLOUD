import sys
sys.path.append('../../../')

from mpl_toolkits.axes_grid1 import make_axes_locatable
import pylab as pl
import subprocess
import numpy as np
import PyECLOUD.mystyle as ms
import os

import scipy.io as sio
import PyECLOUD.myloadmat_to_obj as mlo

pl.close('all')

fontsz = 12

bunch_spacing = 25e-9

first_passage = 1
num_passage = 10

# Please run the "LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_stress_saver" before using this script
folder_sim = '../../testing/tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_em_tracking/'
folder_imgs = 'maps/'
main_outp_filename = 'Pyecltest_angle3D.mat'

passlist = range(first_passage, first_passage + num_passage)

x_beam_pos = 0.
y_beam_pos = 0.

color_beam = 'r'
tbeam_win_length = bunch_spacing * num_passage

########

########

xfield_cut_lim = 7  # 5000.
yfield_cut_lim = 7  # 5000.

denslim = 1.e11

N_dec = 1
lim_dens = [10, 14]
flag_log = True
outp_filname = 'pass%d.avi'%passlist[0]

########

########

qe = 1.60217657e-19
c = 299792458.

i_photog = 0

ms.mystyle_arial(fontsz=fontsz, dist_tick_lab=10)

chamber = mlo.myloadmat_to_obj(folder_sim + '/LHC_chm_ver.mat')
Vx = chamber.Vx
Vy = chamber.Vy

##load output file
obout = mlo.myloadmat_to_obj(folder_sim + '/' + main_outp_filename)
t = obout.t
lam_b1 = obout.lam_t_array
cendens = np.squeeze(obout.cen_density)

#every bunch passage
for pass_ind in passlist:

    filename_rho = 'rho_video/rho_pass%d.mat'%pass_ind
    filename_sc_forces = 'sc_forces/sc_forces_pass%d.mat'%pass_ind
    #new code

    obout_ecl_video = mlo.myloadmat_to_obj(filename_rho)
    obout_sc_forces = mlo.myloadmat_to_obj(filename_sc_forces)

    xg_sc = obout_ecl_video.xg_sc
    yg_sc = obout_ecl_video.yg_sc

    xmin = np.min(xg_sc)
    xmax = np.max(xg_sc)
    ymin = np.min(yg_sc)
    ymax = np.max(yg_sc)

    ix_0 = np.argmin(np.abs(xg_sc - x_beam_pos))
    iy_0 = np.argmin(np.abs(yg_sc - y_beam_pos))

    rho_video = -obout_ecl_video.rho_video / qe
    sc_Fe_norm_video = np.sqrt(np.multiply(obout_sc_forces.Fefx_video,obout_sc_forces.Fefx_video)+np.multiply(obout_sc_forces.Fefy_video,obout_sc_forces.Fefy_video))
    sc_Fes_norm_video = np.sqrt(np.multiply(obout_sc_forces.Fesfx_video,obout_sc_forces.Fesfx_video)+np.multiply(obout_sc_forces.Fesfy_video,obout_sc_forces.Fesfy_video))
    sc_Fb_norm_video = np.sqrt(np.multiply(obout_sc_forces.Fbfx_video,obout_sc_forces.Fbfx_video)+np.multiply(obout_sc_forces.Fbfy_video,obout_sc_forces.Fbfy_video)+np.multiply(obout_sc_forces.Fbfz_video,obout_sc_forces.Fbfz_video))
    t_video_old = np.squeeze(obout_ecl_video.t_video.real)
    t_video = np.squeeze(obout_sc_forces.t_efield_video_new.real)
    b_spac = np.squeeze(obout.b_spac.real)

    (nphotog, _, _) = sc_Fe_norm_video.shape
    (nphotog_big, _, _) = rho_video.shape

    for ii in xrange(0, nphotog, N_dec):
		fig = pl.figure(1, figsize=(9, 9))  # ,figsize=(4.5,6)
		t_curr = t_video[ii]
		ii2, = np.where( t_video_old == t_curr)
		t_curr_old = t_video_old[ii2]
		cendens_curr = np.interp(t_curr, t, cendens)
		lam_b1_curr = np.interp(t_curr, t, lam_b1)
		cendens_curr = np.interp(t_curr, t, cendens)

		print 'Pass %d %d/%d'%(pass_ind, ii, nphotog)

		imm = np.squeeze(rho_video[ii2, :, :])
		imm_Fe_norm = np.squeeze(sc_Fe_norm_video[ii, :, :])
		imm_Fes_norm = np.squeeze(sc_Fes_norm_video[ii, :, :])
		imm_Fb_norm = np.squeeze(sc_Fb_norm_video[ii, :, :])
		immlin = imm.copy()

		if flag_log:

			imm = np.log10((imm))
			imm[immlin < 10**lim_dens[0]] = lim_dens[0]

		pl.subplot(2, 2, 1)
		ax = pl.gca()
		im = ax.imshow(imm.T, cmap=None, norm=None, aspect='auto', interpolation=None,
                 alpha=None, vmin=lim_dens[0], vmax=lim_dens[1], origin='lower', extent=[xmin * 1e3, xmax * 1e3, ymin * 1e3, ymax * 1e3])

		pl.plot(Vx * 1e3, Vy * 1e3, 'y', linewidth=1.5)
		pl.axis('equal')
		pl.xlim(np.min(Vx * 1e3) - np.max(Vx * 1e3) / 10, np.max(Vx * 1e3) + np.max(Vx * 1e3) / 10)
		pl.ylim(np.min(Vy * 1e3) - np.max(Vy * 1e3) / 10, np.max(Vy * 1e3) + np.max(Vy * 1e3) / 10)
		#pl.xlabel('x [mm]')
		#pl.ylabel('y [mm]')

		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.02)
		clb = pl.colorbar(im, cax=cax)
		clb.set_label('log10(e- dens.)')

		pl.subplot(2, 2, 2)
		ax = pl.gca()
		im = ax.imshow(imm_Fe_norm.T, cmap=None, norm=None, aspect='auto', interpolation=None,
                 alpha=None, vmax=0.8*max(map(max, imm_Fe_norm)), origin='lower', extent=[xmin * 1e3, xmax * 1e3, ymin * 1e3, ymax * 1e3])

		pl.plot(Vx * 1e3, Vy * 1e3, 'y', linewidth=1.5)
		pl.axis('equal')
		pl.xlim(np.min(Vx * 1e3) - np.max(Vx * 1e3) / 10, np.max(Vx * 1e3) + np.max(Vx * 1e3) / 10)
		pl.ylim(np.min(Vy * 1e3) - np.max(Vy * 1e3) / 10, np.max(Vy * 1e3) + np.max(Vy * 1e3) / 10)

		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.02)
		clb = pl.colorbar(im, cax=cax)
		clb.set_label('||Fe||   [N]')

		pl.subplot(2, 2, 3)
		ax = pl.gca()
		im = ax.imshow(imm_Fb_norm.T, cmap=None, norm=None, aspect='auto', interpolation=None,
                 alpha=None, vmax=0.8*max(map(max, imm_Fb_norm)), origin='lower', extent=[xmin * 1e3, xmax * 1e3, ymin * 1e3, ymax * 1e3])

		pl.plot(Vx * 1e3, Vy * 1e3, 'y', linewidth=1.5)
		pl.axis('equal')
		pl.xlim(np.min(Vx * 1e3) - np.max(Vx * 1e3) / 10, np.max(Vx * 1e3) + np.max(Vx * 1e3) / 10)
		pl.ylim(np.min(Vy * 1e3) - np.max(Vy * 1e3) / 10, np.max(Vy * 1e3) + np.max(Vy * 1e3) / 10)

		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.02)
		clb = pl.colorbar(im, cax=cax)
		clb.set_label('||Fb||  [N]')

		pl.subplot(2, 2, 4)
		ax = pl.gca()
		im = ax.imshow(imm_Fes_norm.T, cmap=None, norm=None, aspect='auto', interpolation=None,
                 alpha=None, vmax=0.8*max(map(max, imm_Fes_norm)), origin='lower', extent=[xmin * 1e3, xmax * 1e3, ymin * 1e3, ymax * 1e3])

		pl.plot(Vx * 1e3, Vy * 1e3, 'y', linewidth=1.5)
		pl.axis('equal')
		pl.xlim(np.min(Vx * 1e3) - np.max(Vx * 1e3) / 10, np.max(Vx * 1e3) + np.max(Vx * 1e3) / 10)
		pl.ylim(np.min(Vy * 1e3) - np.max(Vy * 1e3) / 10, np.max(Vy * 1e3) + np.max(Vy * 1e3) / 10)

		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.02)
		clb = pl.colorbar(im, cax=cax)
		clb.set_label('||Fe||  [N]  (electrostatic)')

		pl.subplots_adjust(top=0.95, right=0.9, left=0.07, hspace=0.1, wspace=0.7)
		filename = str('maps/Pass%05d_%05d' % (pass_ind, ii)) + '.png'
		pl.savefig(filename, dpi=150)
		pl.clf()

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
           outp_filname)


subprocess.check_call(command)
folderpngname = outp_filname.split('.avi')[0] + '_high_res_pngs'
os.system('mkdir %s'%folderpngname)
os.system('mv *.png ' + folderpngname)

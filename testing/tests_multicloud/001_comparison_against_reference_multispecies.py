import sys, os
BIN = os.path.expanduser("../../../")
sys.path.append(BIN)

import argparse
import pylab as pl
import numpy as np
#from colorsys import hsv_to_rgb
import os
import PyECLOUD.myloadmat_to_obj as mlm
import matplotlib.gridspec as gridspec
import PyECLOUD.mystyle as ms

pl.close('all')

sim_folder = 'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_multispecies'
cloud_name_list = ['electron1', 'electron2', 'electron3']

# check if user provided folder as command line argument
parser = argparse.ArgumentParser()
parser.add_argument('--folder', help='Simulation folder')
parser.add_argument('--angle-dist-func', help='Angular distribution of new MPs relative to surface normal. Introduced in July 2017. Default is 3D.', choices=('2D', '3D'), default='3D')
parser.add_argument('--cloud_list', help='List of cloud names')
args = parser.parse_args()
if args.folder:
    sim_folder = args.folder
if args.cloud_list:
    cloud_name_list = args.cloud_list


# start comparison against reference
ref_folder = '../tests_buildup/LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns/'
curr_folder = sim_folder

folder_plot = sim_folder + '/comparison_plots'


try:
    os.mkdir(folder_plot)
except OSError as err:
    print('Folder not created due to exception:')
    print(err)


myfontsz = 10
titlesz = 10
labelsz = 10
ms.mystyle_arial(fontsz=myfontsz)

color_list = ['k', 'g', 'm']

dict_ref = mlm.myloadmat(ref_folder + '/Pyecltest_angle%s_ref.mat' % args.angle_dist_func) # load dictionary of the reference simulation
dict_curr_list = []
for cloud_name in cloud_name_list:
    dict_curr_list.append(mlm.myloadmat(curr_folder + '/Pyecltest_angle%s_%s.mat' % (args.angle_dist_func, cloud_name)))   # load dictionary of the current simulation


out_var_ref = list(dict_ref.keys())       # returns the list of keys
out_var_curr = list(dict_curr_list[0].keys())

out_var_ref.sort()         # sort the keys in alphabetical order
out_var_curr.sort()

n_pass = 29  # reference passage

var_no_sum_list = ['En_g_hist', 'lam_t_array', 'N_mp_corrected_pass', 'N_mp_impact_pass', 'N_mp_ref_pass', 'sey_test_cos_theta', 'sey_test_del_elast_mat',
                   'sey_test_del_true_mat', 'sey_test_E_impact_eV', 't', 't_En_hist', 't_hist', 't_sc_video', 'xg_hist', 'xg_hist_cos_angle', 'U_sc_eV']

print('Curr Variables are:')
for variab in out_var_curr:
    print(variab, dict_curr_list[0][variab].shape)


for ii, k in enumerate(out_var_curr):
    if '__' in k or k == 'el_dens_at_probes':
        print(('Skipped %s' % k))
        continue

    if len(dict_curr_list[0][k].shape) == 1:  # var is a vector!

        fig = pl.figure(ii)
        fig.patch.set_facecolor('w')
        pl.subplots_adjust(right=0.75)
        pl.title(out_var_curr[ii], fontsize=titlesz)

        #Plot vector for the current simulation
        out_var_curr_tot = dict_curr_list[0][k] * 0
        for i_curr, dict_curr in enumerate(dict_curr_list):
            col = color_list[i_curr]
            pl.plot(dict_curr[k], label=cloud_name_list[i_curr], color=col)
            out_var_curr_tot += dict_curr[k]
        if out_var_curr[ii] not in var_no_sum_list:
            pl.plot(out_var_curr_tot, 'b', label='curr_sim')

        #Plot vector for the reference simulation
        if (k in out_var_ref) and (dict_ref[k].shape != ()):

            pl.plot(dict_ref[k], 'r', label='ref_sim')
            print(ii, k, 'ref_sim')

        else:
            print('%s not  in reference'%k)

        pl.legend(prop={'size': myfontsz}, bbox_to_anchor=(1, 1), loc='best')
        ms.sciy()
        pl.savefig(folder_plot + '/angle%s_%s'%(args.angle_dist_func, k), dpi=300)

    elif len(dict_curr[k].shape) == 2:  # var is a matrix!!!!!!!!!!!!!!!!!

        fig = pl.figure(ii)
        fig.patch.set_facecolor('w')
        pl.subplots_adjust(top=1.2)
        pl.suptitle(out_var_curr[ii], fontsize=titlesz)
        gs1 = gridspec.GridSpec(2, 1)
        gs2 = gridspec.GridSpec(3, 1)

        #Plot matrix for the current simulation
        sp1 = fig.add_subplot(gs1[0])
        sp1.set_title('curr_sim', fontsize=titlesz)
        out_var_curr_tot = dict_curr[k] * 0
        for i_curr, dict_curr in enumerate(dict_curr_list):
            out_var_curr_tot += dict_curr[k]
        pl.pcolormesh(out_var_curr_tot)
        pl.tick_params(labelsize=labelsz)
        cbar = pl.colorbar()
        cbar.ax.tick_params(labelsize=labelsz)
        cbar.formatter.set_powerlimits((0, 0))
        cbar.update_ticks()
        ms.sciy()
        print(ii, k, 'curr_sim')

        try:
            ind_in_mat = n_pass
            if k.startswith('sey_test_'):
                ind_in_mat = 3

            #Plot number of e- for the reference passage
            sp3 = fig.add_subplot(gs2[0])
            out_var_curr_tot = dict_curr[k][ind_in_mat] * 0
            for i_curr, dict_curr in enumerate(dict_curr_list):
                col = color_list[i_curr]
                sp3.plot(dict_curr[k][ind_in_mat], label=cloud_name_list[i_curr], color=col)
                out_var_curr_tot += dict_curr[k][ind_in_mat]
            if out_var_curr[ii] not in var_no_sum_list:
                sp3.plot(out_var_curr_tot, 'b', label='curr_sim')
            sp3.legend(prop={'size': myfontsz}, loc='best')
            sp3.set_title(' num pass equal to [%d]'%ind_in_mat, fontsize=titlesz)
            sp3.tick_params(labelsize=labelsz)
            ms.sciy()

            #Plot number of e- for each slice
            sp4 = fig.add_subplot(gs2[1])
            out_var_curr_tot = dict_curr[k] * 0
            for i_curr, dict_curr in enumerate(dict_curr_list):
                col = color_list[i_curr]
                sp4.plot(np.sum(dict_curr[k], axis=0), label=cloud_name_list[i_curr], color=col)
                out_var_curr_tot += dict_curr[k]
            if out_var_curr[ii] not in var_no_sum_list:
                sp4.plot(np.sum(out_var_curr_tot, axis=0), 'b', label='curr_sim')
            sp4.legend(prop={'size': myfontsz}, loc='best')
            sp4.set_title('e- per slice', fontsize=titlesz)
            sp4.tick_params(labelsize=labelsz)
            ms.sciy()

            #Plot number of e- for each passage
            sp5 = fig.add_subplot(gs2[2])
            out_var_curr_tot = dict_curr[k] * 0
            for i_curr, dict_curr in enumerate(dict_curr_list):
                col = color_list[i_curr]
                sp5.plot(np.sum(dict_curr[k], axis=1), label=cloud_name_list[i_curr], color=col)
                out_var_curr_tot += dict_curr[k]
            if out_var_curr[ii] not in var_no_sum_list:
                sp5.plot(np.sum(out_var_curr_tot, axis=1), 'b', label='curr_sim')
            sp5.legend(prop={'size': myfontsz}, loc='best', )
            sp5.set_title('e- per passage', fontsize=titlesz)
            sp5.tick_params(labelsize=labelsz)
            ms.sciy()

            gs2.tight_layout(fig, rect=[0.45, 0, 1, 1], pad=1.08, h_pad=0.5)

        except IOError as goterror:
            print('Skipped. Got:', goterror)
        except IndexError as goterror:
            print('Skipped. Got:', goterror)

        if (k in out_var_ref) and (dict_ref[k].shape != ()):

            #Plot matrix for the reference simulation
            sp2 = fig.add_subplot(gs1[1])
            sp2.set_title('ref_sim', fontsize=titlesz)
            pl.pcolormesh(dict_ref[k])
            cbar = pl.colorbar()
            cbar.ax.tick_params(labelsize=labelsz)
            cbar.formatter.set_powerlimits((0, 0))
            cbar.update_ticks()
            ms.sciy()
            print(ii, k, 'ref_sim')

            try:

                ind_in_mat = n_pass
                if k.startswith('sey_test_'):
                    ind_in_mat = 3

                #Plot number of e- for the reference passage
                sp3.plot(dict_ref[k][ind_in_mat], 'r', label='ref_sim')
                sp3.legend(prop={'size': myfontsz}, loc='best')
                ms.sciy()

                #Plot number of e- for each slice
                sp4.plot(np.sum(dict_ref[k], axis=0), 'r', label='ref_sim')
                sp4.legend(prop={'size': myfontsz}, loc='best')
                ms.sciy()

                #Plot number of e- for each passage
                sp5.plot(np.sum(dict_ref[k], axis=1), 'r', label='ref_sim')
                sp5.legend(prop={'size': myfontsz}, loc='best')
                ms.sciy()

                gs2.tight_layout(fig, rect=[0.45, 0, 1, 1], pad=1.08, h_pad=1.5)

            except IOError as goterror:
                    print('Skipped. Got:', goterror)
            except IndexError as goterror:
                    print('Skipped. Got:', goterror)

        else:
            print('%s not  in reference'%k)

        gs1.tight_layout(fig, rect=[0, 0, 0.45, 1], pad=1.08)
        top = 0.9
        bottom = max(gs1.bottom, gs2.bottom)

        gs1.update(top=top, bottom=bottom)
        gs2.update(top=top, bottom=bottom)
        pl.savefig(folder_plot + '/angle%s_%s'%(args.angle_dist_func, k), dpi=300)


print('Saved comparison plots in:')
print(folder_plot)

print('In ipython, you may call EOG() to view the results if EOG is installed.')
EOG = lambda : os.system('eog %s/*%s*' % (folder_plot, args.angle_dist_func))
# pl.show()



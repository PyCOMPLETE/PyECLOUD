import sys, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
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

parser = argparse.ArgumentParser()
parser.add_argument('--folder', help='Simulation folder')
args = parser.parse_args()

# Either pass the folder as an argument, or choose it manually here
if args.folder:
    sim_folder = args.folder
else:
    sim_folder = 'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns'
    #sim_folder = 'LHC_ArcDipReal_450GeV_sey1.00_2.5e11ppb_bl_1.00ns_gas_ionization'
    #sim_folder = 'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_change_s_and_E0'
    #sim_folder = 'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_multigrid'
    #sim_folder = 'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns'
    #sim_folder = 'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_circular'
    #sim_folder = 'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_skew'
    #sim_folder = 'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_skew_circular'
    #sim_folder = 'LHC_Sextupole_450GeV_sey1.65_2.5e11ppb_bl_1.00ns'
    #sim_folder = 'LHC_Sextupole_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_skew'



ref_folder = sim_folder
curr_folder = sim_folder

folder_plot = sim_folder+'/comparison_plots'


try:
    os.mkdir(folder_plot)
except OSError as err:
    print 'Folder not created due to exception:'
    print err



myfontsz = 14
ms.mystyle_arial(fontsz=myfontsz)

dict_ref = mlm.myloadmat(ref_folder+'/Pyecltest_ref.mat') # load dictionary of the reference simulation
dict_curr = mlm.myloadmat(curr_folder+'/Pyecltest.mat')  # load dictionary of the current simulation


out_var_ref = dict_ref.keys()       # returns the list of keys
out_var_curr = dict_curr.keys()

out_var_ref.sort()         # sort the keys in alphabetical order
out_var_curr.sort()

n_pass = 35  # reference passage



print'Curr Variables are:'
for variab in out_var_curr:
    print variab, dict_curr[variab].shape


for ii,k in enumerate(out_var_curr):


    if len(dict_curr[k].shape)==1:  # var is a vector!!!!!!!!!!!!!!!


        #Plot vector for the current simulation
        fig=pl.figure(ii)
        pl.subplots_adjust(right=0.75)
        pl.title(out_var_curr[ii])

        pl.plot(dict_curr[k],'b', label='curr_sim')
        print ii,k,'curr_sim'



        if (k in out_var_ref) and (dict_ref[k].shape!=()):


            #Plot vector for the reference simulation
            pl.plot(dict_ref[k],'r', label='ref_sim')
            print ii,k,'ref_sim'


        else:
            print '%s not  in reference'%k



        pl.legend(prop={'size':myfontsz}, bbox_to_anchor=(1, 1),  loc='upper left')
        ms.sciy()
        pl.savefig(folder_plot+'/%s'%k, dpi=300)




    elif len(dict_curr[k].shape)==2:  # var is a matrix!!!!!!!!!!!!!!!!!


        fig=pl.figure(ii)
        pl.subplots_adjust(top=1.2)
        pl.suptitle(out_var_curr[ii])
        gs1 = gridspec.GridSpec(2, 1)
        gs2 = gridspec.GridSpec(3, 1)


        #Plot matrix for the current simulation
        sp1 = fig.add_subplot(gs1[0])
        sp1.set_title('curr_sim')
        pl.pcolormesh(dict_curr[k])
        pl.tick_params(labelsize=10)
        cbar=pl.colorbar()
        cbar.ax.tick_params(labelsize=10)
        cbar.formatter.set_powerlimits((0, 0))
        cbar.update_ticks()
        ms.sciy()
        print ii,k,'curr_sim'

        try:

            #Plot number of e- for the reference passage
            sp3=fig.add_subplot(gs2[0])
            sp3.plot(dict_curr[k][n_pass],'b', label='curr_sim')
            sp3.legend(prop={'size':myfontsz},  loc='upper left')
            sp3.set_title(' num pass equal to [%d]'%n_pass)
            sp3.tick_params(labelsize=10)
            ms.sciy()


            #Plot number of e- for each slice
            sp4=fig.add_subplot(gs2[1])
            sp4.plot(np.sum(dict_curr[k], axis=0),'b', label='curr_sim')
            sp4.legend(prop={'size':myfontsz},  loc='upper left')
            sp4.set_title('e- per slice')
            sp4.tick_params(labelsize=10)
            ms.sciy()


            #Plot number of e- for each passage
            sp5=fig.add_subplot(gs2[2])
            sp5.plot(np.sum(dict_curr[k], axis=1),'b', label='curr_sim')
            sp5.legend(prop={'size':myfontsz},  loc='upper right')
            sp5.set_title('e- per passage')
            sp5.tick_params(labelsize=10)
            ms.sciy()

            gs2.tight_layout(fig,rect=[0.45, 0, 1, 1],pad=1.08, h_pad=0.5)


        except IOError as goterror:
                    print 'Skipped. Got:',  goterror


        if (k in out_var_ref) and  (dict_ref[k].shape!=()):


            #Plot matrix for the reference simulation
            sp2= fig.add_subplot(gs1[1])
            sp2.set_title('ref_sim')
            pl.pcolormesh(dict_ref[k])
            cbar=pl.colorbar()
            cbar.ax.tick_params(labelsize=10)
            cbar.formatter.set_powerlimits((0, 0))
            cbar.update_ticks()
            ms.sciy()
            print ii,k,'ref_sim'

            try:

                #Plot number of e- for the reference passage
                sp3.plot(dict_ref[k][n_pass],'r', label='ref_sim')
                sp3.legend(prop={'size':myfontsz},  loc='upper left')
                ms.sciy()


                #Plot number of e- for each slice
                sp4.plot(np.sum(dict_ref[k], axis=0),'r', label='ref_sim')
                sp4.legend(prop={'size':myfontsz},  loc='upper left')
                ms.sciy()

                #Plot number of e- for each passage
                sp5.plot(np.sum(dict_ref[k], axis=1),'r', label='ref_sim')
                sp5.legend(prop={'size':myfontsz},  loc='upper right')
                ms.sciy()

                gs2.tight_layout(fig,rect=[0.45, 0, 1, 1], pad=1.08,h_pad=1.5)

            except IOError as goterror:
                    print 'Skipped. Got:',  goterror

        else:
            print '%s not  in reference'%k




        gs1.tight_layout(fig,rect=[0, 0, 0.45, 1],pad=1.08)
        top = 0.9
        bottom = max(gs1.bottom, gs2.bottom)

        gs1.update(top=top, bottom=bottom)
        gs2.update(top=top, bottom=bottom)
        pl.savefig(folder_plot+'/%s'%k, dpi=300)

print 'Saved comparison plots in:'
print folder_plot

print 'In ipython, you may call EOG() to view the results if eye of gnome is installed.'
EOG = lambda : os.system('eog %s' % folder_plot)
        #~ #pl.show()



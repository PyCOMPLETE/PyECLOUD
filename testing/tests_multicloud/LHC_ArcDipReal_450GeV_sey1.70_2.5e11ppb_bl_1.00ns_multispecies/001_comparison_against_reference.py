import sys, os
BIN = os.path.expanduser("../../../../")
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

#sim_folder = 'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns'
#sim_folder = 'LHC_ArcDipReal_450GeV_sey1.00_2.5e11ppb_bl_1.00ns_gas_ionization'
#sim_folder = 'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_change_s_and_E0'
#sim_folder = 'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_multigrid'
#sim_folder = 'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns'
#sim_folder = 'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_circular'
#sim_folder = 'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_skew'
#sim_folder = 'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_skew_circular'
#sim_folder = 'LHC_Sextupole_450GeV_sey1.65_2.5e11ppb_bl_1.00ns'
#sim_folder = 'LHC_Sextupole_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_skew'
#sim_folder = 'LHC_ArcDipReal_6500GeV_sey_1.70_1.1e11ppb_b1_1.00ns'
#sim_folder = 'LHC_Drift_6500GeV_sey_1.70_1.1e11ppb_b1_1.00ns'
#sim_folder = 'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_seyfromfile'
#sim_folder = 'LHC_Triplet_Quadrupole_multiple_beams'
#sim_folder = 'CLIC_DRe-_Drift_0.5ns_4.0e9ppb_gas_ionization_ecloud_sey2.0'
#sim_folder = 'CLIC_DRe+_Drift_0.5ns_4.0e9ppb_gas_ionization_ecloud_sey2.0'
#sim_folder = 'CLIC_DRe-_Drift_0.5ns_4.0e9ppb_gas_ionization_ions_A18'
#sim_folder = 'CLIC_DRe+_Drift_0.5ns_4.0e9ppb_gas_ionization_ions_A18'
#sim_folder = 'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns_multispecies'
sim_folder = './'

# check if user provided folder as command line argument
parser = argparse.ArgumentParser()
parser.add_argument('--folder', help='Simulation folder')
parser.add_argument('--angle-dist-func', help='Angular distribution of new MPs relative to surface normal. Introduced in July 2017. Default is 3D.', choices=('2D', '3D'), default='3D')
args = parser.parse_args()
if args.folder:
    sim_folder = args.folder


# start comparison against reference
ref_folder = '../LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns/'
curr_folder = sim_folder

folder_plot = sim_folder+'/comparison_plots'


try:
    os.mkdir(folder_plot)
except OSError as err:
    print 'Folder not created due to exception:'
    print err



myfontsz = 10
titlesz = 10
labelsz = 10
ms.mystyle_arial(fontsz=myfontsz)

dict_ref = mlm.myloadmat(ref_folder+'/Pyecltest_angle%s_ref.mat' % args.angle_dist_func) # load dictionary of the reference simulation
dict_curr = mlm.myloadmat(curr_folder+'/Pyecltest_angle%s_electron1.mat' % args.angle_dist_func)   # load dictionary of the current simulation
dict_curr2 = mlm.myloadmat(curr_folder+'/Pyecltest_angle%s_electron2.mat' % args.angle_dist_func)   # load dictionary of the current simulation
dict_curr3 = mlm.myloadmat(curr_folder+'/Pyecltest_angle%s_electron3.mat' % args.angle_dist_func)   # load dictionary of the current simulation


out_var_ref = dict_ref.keys()       # returns the list of keys
out_var_curr = dict_curr.keys()

out_var_ref.sort()         # sort the keys in alphabetical order
out_var_curr.sort()

n_pass = 35  # reference passage



print'Curr Variables are:'
for variab in out_var_curr:
    print variab, dict_curr[variab].shape


for ii,k in enumerate(out_var_curr):
    if '__' in k or k == 'el_dens_at_probes':
        print('Skipped %s' % k)
        continue


    if len(dict_curr[k].shape)==1:  # var is a vector!


        #Plot vector for the current simulation
        fig=pl.figure(ii)
        fig.patch.set_facecolor('w')
        pl.subplots_adjust(right=0.75)
        pl.title(out_var_curr[ii], fontsize=titlesz)

        pl.plot(dict_curr[k],'k', label='cloud 1')
        pl.plot(dict_curr2[k],'g', label='cloud 2')
        pl.plot(dict_curr3[k],'m', label='cloud 3')
        pl.plot(dict_curr[k] + dict_curr2[k] + dict_curr3[k],'b', label='curr_sim')
        print ii,k,'curr_sim'



        if (k in out_var_ref) and (dict_ref[k].shape!=()):


            #Plot vector for the reference simulation
            pl.plot(dict_ref[k],'r', label='ref_sim')
            print ii,k,'ref_sim'


        else:
            print '%s not  in reference'%k



        pl.legend(prop={'size':myfontsz}, bbox_to_anchor=(1, 1),  loc=0)
        ms.sciy()
        pl.savefig(folder_plot+'/angle%s_%s'%(args.angle_dist_func, k), dpi=300)




    elif len(dict_curr[k].shape)==2:  # var is a matrix!!!!!!!!!!!!!!!!!


        fig=pl.figure(ii)
        fig.patch.set_facecolor('w')
        pl.subplots_adjust(top=1.2)
        pl.suptitle(out_var_curr[ii], fontsize=titlesz)
        gs1 = gridspec.GridSpec(2, 1)
        gs2 = gridspec.GridSpec(3, 1)


        #Plot matrix for the current simulation
        sp1 = fig.add_subplot(gs1[0])
        sp1.set_title('curr_sim', fontsize=titlesz)
        pl.pcolormesh(dict_curr[k] + dict_curr2[k] + dict_curr3[k])
        pl.tick_params(labelsize=labelsz)
        cbar=pl.colorbar()
        cbar.ax.tick_params(labelsize=labelsz)
        cbar.formatter.set_powerlimits((0, 0))
        cbar.update_ticks()
        ms.sciy()
        print ii,k,'curr_sim'

        try:

            #Plot number of e- for the reference passage
            sp3=fig.add_subplot(gs2[0])
            sp3.plot(dict_curr[k][n_pass] + dict_curr2[k][n_pass] + dict_curr3[k][n_pass],'b')
            sp3.plot(dict_curr[k][n_pass],'k')
            sp3.plot(dict_curr2[k][n_pass],'g')
            sp3.plot(dict_curr3[k][n_pass],'m')
            sp3.legend(prop={'size':myfontsz},  loc=0)
            sp3.set_title(' num pass equal to [%d]'%n_pass, fontsize=titlesz)
            sp3.tick_params(labelsize=labelsz)
            ms.sciy()


            #Plot number of e- for each slice
            sp4=fig.add_subplot(gs2[1])
            sp4.plot(np.sum(dict_curr[k], axis=0) + np.sum(dict_curr2[k], axis=0) + np.sum(dict_curr3[k], axis=0),'b', label='curr_sim')
            sp4.plot(np.sum(dict_curr[k], axis=0),'k', label='cloud1')
            sp4.plot(np.sum(dict_curr2[k], axis=0),'g', label='cloud2')
            sp4.plot(np.sum(dict_curr3[k], axis=0),'m', label='cloud3')
            sp4.legend(prop={'size':myfontsz},  loc=0, frameon=False)
            sp4.set_title('e- per slice', fontsize=titlesz)
            sp4.tick_params(labelsize=labelsz)
            ms.sciy()


            #Plot number of e- for each passage
            sp5=fig.add_subplot(gs2[2])
            sp5.plot(np.sum(dict_curr[k], axis=1) + np.sum(dict_curr2[k], axis=1) + np.sum(dict_curr3[k], axis=1),'b')
            sp5.plot(np.sum(dict_curr[k], axis=1),'k')
            sp5.plot(np.sum(dict_curr2[k], axis=1),'g')
            sp5.plot(np.sum(dict_curr3[k], axis=1),'m')
            sp5.legend(prop={'size':myfontsz},  loc=0, )
            sp5.set_title('e- per passage', fontsize=titlesz)
            sp5.tick_params(labelsize=labelsz)
            ms.sciy()

            gs2.tight_layout(fig,rect=[0.45, 0, 1, 1],pad=1.08, h_pad=0.5)


        except IOError as goterror:
                    print 'Skipped. Got:',  goterror


        if (k in out_var_ref) and  (dict_ref[k].shape!=()):


            #Plot matrix for the reference simulation
            sp2= fig.add_subplot(gs1[1])
            sp2.set_title('ref_sim', fontsize=titlesz)
            pl.pcolormesh(dict_ref[k])
            cbar=pl.colorbar()
            cbar.ax.tick_params(labelsize=labelsz)
            cbar.formatter.set_powerlimits((0, 0))
            cbar.update_ticks()
            ms.sciy()
            print ii,k,'ref_sim'

            try:

                #Plot number of e- for the reference passage
                sp3.plot(dict_ref[k][n_pass],'r')
                sp3.legend(prop={'size':myfontsz},  loc=0)
                ms.sciy()


                #Plot number of e- for each slice
                sp4.plot(np.sum(dict_ref[k], axis=0),'r', label='ref_sim')
                sp4.legend(prop={'size':myfontsz},  loc=0, frameon=False)
                ms.sciy()

                #Plot number of e- for each passage
                sp5.plot(np.sum(dict_ref[k], axis=1),'r')
                sp5.legend(prop={'size':myfontsz},  loc=0)
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
        pl.savefig(folder_plot+'/angle%s_%s'%(args.angle_dist_func, k), dpi=300)

print 'Saved comparison plots in:'
print folder_plot

print 'In ipython, you may call EOG() to view the results if EOG is installed.'
EOG = lambda : os.system('eog %s/*%s*' % (folder_plot, args.angle_dist_func))
# pl.show()


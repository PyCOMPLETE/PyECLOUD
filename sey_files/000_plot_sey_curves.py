from __future__ import division
import os
import re
import numpy as np
import matplotlib.pyplot as plt

import LHCMeasurementTools.mystyle as ms

import PyECLOUD.sec_emission_model_from_file as sem
import PyECLOUD.sec_emission_model_ECLOUD as sem_e

plt.close('all')
ms.mystyle(12)

plt.rcParams['lines.markersize'] = 3


regex_label = re.compile('sample_(.*cm_.*)\.txt$')
main_dir = os.path.abspath(os.path.dirname(__file__)) + '/SEY-LE_SEY'
all_files = os.listdir(main_dir)
plot_files = sorted(filter(lambda x: 'merged' in x and x.endswith('.txt'), all_files))

xx_lin = np.exp(np.linspace(np.log(0.1), np.log(1.8e3), int(1e3)))

sey_model_e = sem_e.SEY_model_ECLOUD(330, 1.85, 0.7)

fig = ms.figure('SEY from files')

sp1 = plt.subplot(2,2,1)
sp2 = plt.subplot(2,2,2)
sp3 = plt.subplot(2,2,3)
sp4 = plt.subplot(2,2,4)
sp1.set_title('Low energy')
sp2.set_title('High energy')
sp3.set_title('Full range (Logarithmic)')
sp4.set_title('Implementation (Logarithmic)')
for sp in sp1, sp2:
    sp.set_xlabel('Energy [eV]')
    sp.set_ylabel('SEY')

for file_ in plot_files:
    label = regex_label.search(file_).groups()[0]
    label_full = label

    if '3.5' in file_:
        color = 'b'
    else:
        color = 'g'

    sef = sem.SEY_model_from_file(main_dir + '/' + file_, range_extrapolate_right=300, delta_e=0.1, flag_factor_costheta=False, max_sey=False)

    xx, yy = sef.energy_eV_0, sef.sey_parameter_0
    sp2.plot(xx, yy, '.', label=label, ls='None', color=color)
    sp3.semilogx(xx, yy,'.', label=label_full, ls='None', color=color)
    sp1.plot(xx, yy,'.', label=label_full, ls='None', color=color)

    xx_lin_r = xx_lin.copy()
    np.random.shuffle(xx_lin_r)
    yy_r = sef.SEY_process(1, xx_lin_r, 1, None)[0]
    sp4.semilogx(xx_lin_r, yy_r,'.', label=label_full, ls='None', color=color)


yy = sey_model_e.SEY_process(1, xx_lin, 1, None)[0]

for sp_ in sp2, sp3, sp4:
    sp_.plot(xx_lin, yy, label='ECLOUD model', color='r')


xx = np.linspace(1, 80, int(1e3))
yy = sey_model_e.SEY_process(1, xx, 1, None)[0]
sp1.plot(xx, yy, label='ECLOUD model', color='r')

sp1.set_xlim(0, 80)

for sp in sp1, sp2, sp3, sp4:
    sp.legend(loc='upper left', bbox_to_anchor=(1,1))



sey_vect = [1.1, 1.4, 1.7, False][::-1]
ms.figure('SEY from files, max SEY')

for file_ctr, file_ in enumerate(plot_files):
    sp_ctr = file_ctr + 1
    sp = plt.subplot(2,2,sp_ctr)
    sp.grid(True)
    sp.set_xlabel('Energy [eV]')
    sp.set_ylabel('SEY')
    sp.set_title(file_)

    label = regex_label.search(file_).groups()[0]
    label_full = label

    for sey_ctr, max_sey in enumerate(sey_vect):
        color = ms.colorprog(sey_ctr, sey_vect)


        sef = sem.SEY_model_from_file(main_dir + '/' + file_, range_extrapolate_right=300, delta_e=0.1, flag_factor_costheta=False, max_sey=max_sey)

        if max_sey is False:
            max_sey = sef.sey_parameter_0.max()

        sey_model_e = sem_e.SEY_model_ECLOUD(330, max_sey, 0.7)

        xx, yy = sef.energy_eV_0, sef.sey_parameter_0


        sp.semilogx(xx, yy, marker='o',label='%.2f' % max_sey, color=color)

        yy2 = sey_model_e.SEY_process(1, xx, 1, None)[0]
        sp.semilogx(xx, yy2,ls='--', color=color)

    sp.legend(title='Max SEY')


plt.show()


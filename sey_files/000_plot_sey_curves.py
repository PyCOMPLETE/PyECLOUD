import os
import re
import numpy as np
import matplotlib.pyplot as plt

import LHCMeasurementTools.mystyle as ms

import PyECLOUD.sec_emission_model_from_file as sem
import PyECLOUD.sec_emission_model_ECLOUD as sem_e

plt.close('all')
ms.mystyle(12)


regex_label = re.compile('sample_(.*cm_.*)\.txt$')
main_dir = './SEY-LE_SEY'
all_files = os.listdir(main_dir)
plot_files = sorted(filter(lambda x: 'V2' in x and x.endswith('.txt'), all_files))

xx_lin = np.linspace(1, 1.8e3, 1e3)

sey_model_e = sem_e.SEY_model_ECLOUD(330, 1.85, 0.7)

fig = ms.figure('SEY from files')

sp1 = plt.subplot(2,2,1)
sp2 = plt.subplot(2,2,2)
sp3 = plt.subplot(2,2,3)
sp4 = plt.subplot(2,2,4)
sp1.set_title('Low energy')
sp2.set_title('High energy')
sp3.set_title('Full range')
sp4.set_title('Implementation')
for sp in sp1, sp2:
    sp.set_xlabel('Energy [eV]')
    sp.set_ylabel('SEY')

for file_ in plot_files:
    print file_
    label = regex_label.search(file_).groups()[0]

    if 'LE' in file_:
        sp = sp1
        label_full = 'LE ' + label
        ls = '--'
    else:
        sp = sp2
        label_full = label
        ls = '-'

    if '3.5' in file_:
        color = 'b'
    else:
        color = 'g'

    sef = sem.SEY_from_file(main_dir + '/' + file_, 0.7)

    print '%.2f eV' % sef.work_function

    xx, yy = sef.energy_eV, sef.sey_parameter
    sp.plot(xx, yy, label=label, ls=ls, color=color)
    sp3.semilogx(xx, yy, label=label_full, ls=ls, color=color)

    yy = sef.SEY_process(1, xx_lin, 1, None)[0]
    sp4.semilogx(xx_lin, yy, label=label_full, ls=ls, color=color)


yy = sey_model_e.SEY_process(1, xx_lin, 1, None)[0]

sp2.plot(xx_lin, yy, label='ECLOUD model', color='r')
sp3.plot(xx_lin, yy, label='ECLOUD model', color='r')


xx = np.linspace(0, 80, 1e3)
yy = sey_model_e.SEY_process(1, xx, 1, None)[0]
sp1.plot(xx, yy, label='ECLOUD model', color='r')

for sp in sp1, sp2, sp3, sp4:
    sp.legend(loc='upper left', bbox_to_anchor=(1,1))

plt.show()


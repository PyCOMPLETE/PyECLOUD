import sys
import os
BIN = os.path.expanduser("../../../")
sys.path.append(BIN)

import numpy as np
import matplotlib.pyplot as pl
import seaborn as sns

import PyECLOUD.myloadmat_to_obj as mlm

# Load files
sim_folder = 'LHC_ArcDriftReal_450GeV_electron_n2_sey1.75_P3.125e+06/'
file_name_sigma = 'cross_section_n2_15.59_20000.00_N1000_log'

sigma_ref = mlm.myloadmat(sim_folder + file_name_sigma)
sigma_extr = mlm.myloadmat(sim_folder + file_name_sigma + '_extracted')

# Define plot environment
sns.set_context('poster', rc={'lines.linewidth': 2.5})
sns.set_style('whitegrid', {'grid.linestyle': '--', 'legend.frameon': True})

pl.close('all')

pl.figure(figsize=(8,6))
pl.semilogx(sigma_extr['energy_eV'], sigma_extr['sigma_cm2_interp'], '.r', label='interpolated')
pl.semilogx(sigma_ref['energy_eV'], sigma_ref['cross_section_cm2'], '--k', label='input')

pl.legend()
pl.xlabel('Electron energy [eV]')
pl.ylabel('Ionization cross-section [cm$^2$]')
pl.subplots_adjust(right=0.93, bottom=0.14)
pl.savefig('Cross_section_interp.png', dpi=200)

pl.figure(figsize=(8,6))
pl.semilogx(sigma_extr['energy_eV'], sigma_extr['sigma_cm2_sampled_electron'], '.', label='sampled with e$^-$')
pl.semilogx(sigma_extr['energy_eV'], sigma_extr['sigma_cm2_sampled_n2'], '.', label='sampled with N$_2\!$$^+$')

pl.semilogx(sigma_ref['energy_eV'], sigma_ref['cross_section_cm2'], '--k', label='input')

pl.legend()
pl.xlabel('Electron energy [eV]')
pl.ylabel('Ionization cross-section [cm$^2$]')
pl.subplots_adjust(right=0.93, bottom=0.14)
pl.savefig('Cross_section_sampled.png', dpi=200)

pl.figure(figsize=(5,7))
pl.semilogy(sigma_ref['cross_section_cm2']*1e16, sigma_ref['energy_eV'], 'k')
pl.ylabel('Electron energy [eV]')
pl.xlabel('Cross-section [10$^{-16}$ cm$^2$]')
pl.ylim(1, 1e4)
pl.subplots_adjust(left=0.20, bottom=0.12)
pl.savefig('Cross_section_vs_energy.png', dpi=200)


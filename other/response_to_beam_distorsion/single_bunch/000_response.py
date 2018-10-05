import sys
sys.path.append('../../../../')

from scipy.constants import c as clight, e as qe
import matplotlib.pyplot as plt
import numpy as np

import PyECLOUD.PyEC4PyHT as PyEC4PyHT
from PyHEADTAIL.particles.slicing import UniformBinSlicer
import PyECLOUD.mystyle as ms

from LHC_custom import LHC


machine_configuration = 'HLLHC-injection'
machine = LHC(n_segments = 1, machine_configuration = machine_configuration)

bunch = machine.generate_6D_Gaussian_bunch(n_macroparticles=300000,
                intensity=1.15e11, epsn_x=2.5e-6, epsn_y=2.5e-6, sigma_z=0.11)
                
bunch.x[bunch.z<5e-2] += 1e-3
                
                

ecloud_ele = PyEC4PyHT.Ecloud(slice_by_slice_mode=True,
            L_ecloud=1., slicer=None, 
            Dt_ref=25e-12, pyecl_input_folder='pyecloud_config',
       ) 


n_slices = 150
z_cut = 2.5e-9/2*clight

slicer = UniformBinSlicer(n_slices = n_slices, z_cuts=(-z_cut, z_cut))
slices_list_for_map = bunch.extract_slices(slicer)
        
ecloud_ele.save_ele_distributions_last_track = True
ecloud_ele.save_ele_field = True
ecloud_ele._reinitialize()


z_centers = []
# sim bunch-ecloud interaction
for ss in slices_list_for_map[::-1]:
    z_centers.append(ss.slice_info['z_bin_center'])
    ecloud_ele.track(ss)
    
ecloud_ele._finalize()


plt.close('all')
ms.mystyle_arial(fontsz=14, dist_tick_lab=5)

plt.figure(4)
y_beam_offset = 0.
vmax = 3e12
i_y = np.argmin(np.abs(ecloud_ele.spacech_ele.yg-y_beam_offset))
plt.pcolormesh(np.array(z_centers), ecloud_ele.spacech_ele.xg, -1/qe*ecloud_ele.rho_ele_last_track[:,:, i_y].T, vmax=vmax,
    shading='Gouraud'
    )
plt.ylim(-4e-3, 4e-3)
plt.colorbar()
plt.show()

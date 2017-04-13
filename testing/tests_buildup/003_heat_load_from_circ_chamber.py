import scipy.io as sio
import numpy as np

sim_folders = [
    'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_skew_circular',
    'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns_circular',
]

hls = []
for sim in sim_folders:
    mat = sio.loadmat(sim + '/Pyecltest.mat')
    hl = np.sum(mat['energ_eV_impact_hist'])
    hls.append(hl)
    print hl

print ((hls[1] - hls[0]) / hls[0])


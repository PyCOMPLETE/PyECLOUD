import numpy as np
import NAFFlib
import sys


def tune_analysis(x_i, xp_i, y_i, yp_i):
    n_turns = x_i.shape[1]
    macroparticlenumber = x_i.shape[0]

    qx_i = np.empty(macroparticlenumber)
    qy_i = np.empty(macroparticlenumber)

    print('analysing particle spectra ... this may take some time.')
    for p_idx in range(macroparticlenumber):
        
        qx_i[p_idx] = NAFFlib.get_tune(x_i[p_idx, :])
        qy_i[p_idx] = NAFFlib.get_tune(y_i[p_idx, :])

        sys.stdout.write('\rparticle %d'%p_idx)

    x_centroid = np.mean(x_i, axis=0)
    y_centroid = np.mean(y_i, axis=0)

    #print x_centroid.shape
    qx_centroid = NAFFlib.get_tune(x_centroid) 
    qy_centroid = NAFFlib.get_tune(y_centroid) 

    return qx_i, qy_i, qx_centroid, qy_centroid



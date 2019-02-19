import numpy as np
import pylab as pl
import scipy.io as sio

from PyECLOUD import geom_impact_rect_fast_impact as girfi
from PyECLOUD import geom_impact_poly_fast_impact as gipfi

n_theta = 90

x_ellip = 2.3e-2 
y_ellip = 2.3e-2
x_rect = 3e-2
y_rect = 1.8e-2

theta_unif = np.linspace(0, 2*np.pi, n_theta+1)[:-1]

Vx_unif = x_ellip*np.cos(theta_unif);
Vy_unif = y_ellip*np.sin(theta_unif);

Vx_unif[Vx_unif>x_rect]=x_rect;
Vx_unif[Vx_unif<-x_rect]=-x_rect;

Vy_unif[Vy_unif>y_rect]=y_rect;
Vy_unif[Vy_unif<-y_rect]=-y_rect;

x_sem_ellip_insc = 0.98*np.min([x_ellip , x_rect]);
y_sem_ellip_insc = 0.98*np.min([y_ellip , y_rect]);

mask_patch = np.logical_and(Vx_unif>=5.e-3, Vx_unif <= 1e-2)
mask_patch = np.logical_and(Vy_unif>0, mask_patch)

# Regular chamber surface
del_max_segments = Vx_unif*0. +1.35
Emax_segments = Vx_unif*0. + 332.
R0_segments = Vx_unif*0. + 0.7,

flag_charging = Vx_unif*0. 
Q_max_segments = Vx_unif*0. 
EQ_segments = Vx_unif*0. 

# Charging patch
del_max_segments[mask_patch] = 1.9 
flag_charging[mask_patch] = 1. 
Q_max_segments[mask_patch] = 0.1*1e-12*1e6 # 1e-12 C/mm^2 
EQ_segments[mask_patch] = 0.0000000001#20. 


sio.savemat('chamber.mat',{\
'Vx':Vx_unif, \
'Vy':Vy_unif, \
'x_sem_ellip_insc':x_sem_ellip_insc, 
'y_sem_ellip_insc':y_sem_ellip_insc,
'x_aper':np.max(Vx_unif),
'y_aper':np.max(Vy_unif),

'del_max_segments': del_max_segments,
'R0_segments': R0_segments,
'Emax_segments': Emax_segments,

'flag_charging': flag_charging,
'Q_max_segments': Q_max_segments,
'EQ_segments': EQ_segments,

}, oned_as='row')
        
pl.close('all')   
        

pl.figure(2)
pl.plot(Vx_unif, Vy_unif, 'b.-')
pl.plot(Vx_unif[mask_patch], Vy_unif[mask_patch], 'r.-')

pl.plot(x_sem_ellip_insc*np.cos(theta_unif), y_sem_ellip_insc*np.sin(theta_unif), 'g')
pl.axis('equal')

pl.show()
    

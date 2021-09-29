import numpy as np
import pylab as pl
import scipy.io as sio

from PyECLOUD import geom_impact_rect_fast_impact as girfi
from PyECLOUD import geom_impact_poly_fast_impact as gipfi

n_theta = 150

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

Vx_rep = np.concatenate((Vx_unif, [Vx_unif[0]]))
Vx_mid = 0.5*(Vx_rep[:-1]+Vx_rep[1:])

sio.savemat('chamber.mat',{\
'Vx':Vx_unif, \
'Vy':Vy_unif, \
'x_sem_ellip_insc':x_sem_ellip_insc, 
'y_sem_ellip_insc':y_sem_ellip_insc,
'x_aper':np.max(Vx_unif),
'y_aper':np.max(Vy_unif),
}, oned_as='row')
        
pl.close('all')   
        

pl.figure(2)
pl.plot(Vx_unif, Vy_unif, 'b.-')

pl.plot(x_sem_ellip_insc*np.cos(theta_unif), y_sem_ellip_insc*np.sin(theta_unif), 'g')
pl.axis('equal')

pl.show()
    

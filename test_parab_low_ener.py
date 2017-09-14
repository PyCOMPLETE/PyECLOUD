from . import sec_emission_model_cos_low_ener as sem_cle
from . import sec_emission_model_flat_low_ener as sem_fle
from . import sec_emission_model_ECLOUD as sem_ECLOUD
import numpy as np
import pylab as pl
E_vect=np.linspace(0.,1000., 20000)

delta_vec_cle, ref_frac_vec_le=sem_cle.yield_fun2(E_vect,1.,332.,1.,0.7)
delta_vec_fle, ref_frac_vec_fle=sem_fle.yield_fun2(E_vect,1.,332.,1.,0.7)
delta_vec_ecl, ref_frac_vec_ecl=sem_ECLOUD.yield_fun2(E_vect,1.,332.,1.,0.7)

pl.close('all')
pl.plot(E_vect, delta_vec_cle)
pl.plot(E_vect, delta_vec_fle)
pl.plot(E_vect, delta_vec_ecl)
pl.show()

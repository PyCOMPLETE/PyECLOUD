#----------------------------------------------------------------------
#                                                                      
#                           CERN                                       
#                                                                      
#     European Organization for Nuclear Research                       
#                                                                      
#     
#     This file is part of the code:
#                                                                                           
#                          PyECLOUD Version 4.32                   
#                  
#                                                                       
#     Author and contact:   Giovanni IADAROLA 
#                           BE-ABP Group                               
#                           CERN                                       
#                           CH-1211 GENEVA 23                          
#                           SWITZERLAND  
#                           giovanni.iadarola@cern.ch                  
#                                                                      
#                contact:   Giovanni RUMOLO                            
#                           BE-ABP Group                               
#                           CERN                                      
#                           CH-1211 GENEVA 23                          
#                           SWITZERLAND  
#                           giovanni.rumolo@cern.ch                    
#                                                                      
#
#                                                                      
#     Copyright  CERN,  Geneva  2011  -  Copyright  and  any   other   
#     appropriate  legal  protection  of  this  computer program and   
#     associated documentation reserved  in  all  countries  of  the   
#     world.                                                           
#                                                                      
#     Organizations collaborating with CERN may receive this program   
#     and documentation freely and without charge.                     
#                                                                      
#     CERN undertakes no obligation  for  the  maintenance  of  this   
#     program,  nor responsibility for its correctness,  and accepts   
#     no liability whatsoever resulting from its use.                  
#                                                                      
#     Program  and documentation are provided solely for the use  of   
#     the organization to which they are distributed.                  
#                                                                      
#     This program  may  not  be  copied  or  otherwise  distributed   
#     without  permission. This message must be retained on this and   
#     any other authorized copies.                                     
#                                                                      
#     The material cannot be sold. CERN should be  given  credit  in   
#     all references.                                                  
#----------------------------------------------------------------------






from numpy import array, cross, sum, squeeze
import scipy.io as sio
from boris_cython import boris_step_multipole



  
class pusher_Boris_multipole():
    
    def __init__(self, Dt, N_sub_steps=1, B_multip = []):
        
        print "Tracker: Boris multipole"
        
        self.N_sub_steps = N_sub_steps
        self.Dt = Dt
        
        if Dt is None or N_sub_steps is None:
             self.Dtt = None
        else:
            self.Dtt = Dt / float(N_sub_steps)
        
        self.B_multip = array(B_multip)
        
        if len(self.B_multip)==0:
			self.B_multip=array([0])
        
        if len(self.B_multip)>2:
			raise ValueError('Higher multipoles not yet implemented!')
                   
        print "N_subst_init=%d"% self.N_sub_steps
        
    #@profile
    def step(self, MP_e, Ex_n,Ey_n, Ez_n=0.):
        
        
        if MP_e.N_mp>0:
            
            xn1 = MP_e.x_mp[0:MP_e.N_mp]
            yn1 = MP_e.y_mp[0:MP_e.N_mp]
            zn1 = MP_e.z_mp[0:MP_e.N_mp]
            vxn1 = MP_e.vx_mp[0:MP_e.N_mp]
            vyn1 = MP_e.vy_mp[0:MP_e.N_mp]
            vzn1 = MP_e.vz_mp[0:MP_e.N_mp]
        
            if  Ez_n!=0.:
                raise ValueError('Oooops! Not implemented....')
                

            boris_step_multipole(self.N_sub_steps, self.Dtt, self.B_multip, 
					  xn1, yn1, zn1, vxn1, vyn1, vzn1,
					  Ex_n, Ey_n)

        return MP_e
        
    def stepcustomDt(self, MP_e, Ex_n,Ey_n, Ez_n=0., Dt_substep=None, N_sub_steps=None):
        
        if MP_e.N_mp>0:
            
            xn1 = MP_e.x_mp[0:MP_e.N_mp]
            yn1 = MP_e.y_mp[0:MP_e.N_mp]
            zn1 = MP_e.z_mp[0:MP_e.N_mp]
            vxn1 = MP_e.vx_mp[0:MP_e.N_mp]
            vyn1 = MP_e.vy_mp[0:MP_e.N_mp]
            vzn1 = MP_e.vz_mp[0:MP_e.N_mp]
            
            
            if  Ez_n!=0.:
                raise ValueError('Oooops! Not implemented....')
                

            boris_step_multipole(N_sub_steps, Dt_substep, self.B_multip, 
					  xn1, yn1, zn1, vxn1, vyn1, vzn1,
					  Ex_n, Ey_n)

        return MP_e       
     

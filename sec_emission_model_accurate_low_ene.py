#----------------------------------------------------------------------
#                                                                      
#                           CERN                                       
#                                                                      
#     European Organization for Nuclear Research                       
#                                                                      
#     
#     This file is part of the code:
#                                                                      		            
#		           PyECLOUD Version 5.5.2                     
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

from numpy import sqrt, exp
from numpy.random import rand

def yield_fun3(E,costheta,Emax, s, del_max,R0,E0):

    del_max_tilde=del_max*exp(0.5*(1-costheta))
    E_max_tilde=Emax*(1+0.7*(1-costheta))
    
    x=E/E_max_tilde
    
    del_true_sec=del_max_tilde*(s*x)/(s-1+x**s)
    del_reflected=R0*((sqrt(E)-sqrt(E+E0))/(sqrt(E)+sqrt(E+E0)))**2
    
    delta=del_true_sec+del_reflected
    
    ref_frac=del_reflected/delta
    
    
    return del_true_sec, del_reflected, ref_frac



class SEY_model_acc_low_ene:
    def __init__(self, Emax, s_model, del_max, R0, E0=150):
            self.Emax = Emax
            self.s_model = s_model
            self.del_max = del_max
            self.R0 = R0
            self.E0 = E0
            
    def SEY_process(self,nel_impact,E_impact_eV, costheta_impact, i_impact):
            
            del_true_sec, del_reflected, ref_frac = yield_fun3(E_impact_eV,costheta_impact,self.Emax, self.s_model, self.del_max,self.R0, E0=self.E0);
            ref_prob=del_reflected;
            beta_ts=del_true_sec/(1.-del_reflected);
            
            flag_elast=(rand(len(ref_prob))<ref_prob);
            flag_truesec=~(flag_elast);
    
            
            nel_emit=nel_impact;
            nel_emit[flag_truesec]=nel_impact[flag_truesec]*beta_ts[flag_truesec]; 
            
            return  nel_emit, flag_elast, flag_truesec   
        

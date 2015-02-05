#----------------------------------------------------------------------
#                                                                      
#                           CERN                                       
#                                                                      
#     European Organization for Nuclear Research                       
#                                                                      
#     
#     This file is part of the code:
#                                                                      		            
#		           PyECLOUD Version 4.14                     
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

from numpy import sqrt, exp, cos,pi
from numpy.random import rand

def yield_fun2(E,costheta,Emax,del_max,R0):
    
    s=1.35;
    E0=30.;
    del_max_tilde=del_max*exp(0.5*(1.-costheta));
    E_max_tilde=Emax*(1.+0.7*(1.-costheta));

    x=E/E_max_tilde;
    
    true_sec=del_max_tilde*(s*x)/(s-1.+x**s);
    reflected=0.*true_sec
    mask_ref=E<E0
    #reflected[mask_ref]=R0*(1.-(E[mask_ref]/E0)**2.);
    reflected[mask_ref]=R0*(cos(0.5*pi*E[mask_ref]/E0)**2.);
    
    delta=true_sec+reflected;
    
    ref_frac=reflected/delta;
       
    return delta, ref_frac


class SEY_model_cos_le:
    def __init__(self, Emax,del_max,R0):
            self.Emax = Emax
            self.del_max = del_max
            self.R0 = R0
            print 'Secondary emission model: Cosine Low Energy '
            
    def SEY_process(self,nel_impact,E_impact_eV, costheta_impact, i_impact):
            yiel, ref_frac=yield_fun2(E_impact_eV,costheta_impact,self.Emax,self.del_max,self.R0);
            flag_elast=(rand(len(ref_frac))<ref_frac);
            flag_truesec=~(flag_elast);
            nel_emit=nel_impact*yiel;
            
            return  nel_emit, flag_elast, flag_truesec   
        

#----------------------------------------------------------------------
#                                                                      
#                           CERN                                       
#                                                                      
#     European Organization for Nuclear Research                       
#                                                                      
#     
#     This file is part of the code:
#                                                                      		            
#		           PyECLOUD Version 4.35                     
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



from numpy import *
from numpy.random import lognormal
from numpy.random import randn

def hilleret_model2(switch_no_increase_energy, Ngen, sigmafit, mufit, E_th, En_impact_eV, thresh_low_energy):
    
    if switch_no_increase_energy==0:
        en_eV=lognormal(mufit,sigmafit,Ngen);
        flag_above_th=(en_eV>E_th);
        Nabove_th=sum(flag_above_th);
        
        while Nabove_th>0:
            en_eV[flag_above_th]=lognormal(mufit,sigmafit,Nabove_th);
            flag_above_th=(en_eV>E_th);
            Nabove_th=sum(flag_above_th);
  
    
    
    
    if switch_no_increase_energy==1:
        
        en_eV=0*En_impact_eV;
        
        
        flag_low_energy = En_impact_eV<thresh_low_energy
        flag_high_energy = ~(flag_low_energy)
        N_low_ene=sum(flag_low_energy)
        N_high_ene=sum(flag_high_energy)
        
        #generate low energy
        en_eV_le=randn(N_low_ene)   #in eV   
        flag_negat=logical_or(en_eV_le<0., en_eV_le>4.)
        N_neg=sum(flag_negat);
        while(N_neg>0):
            en_eV_le[flag_negat]=randn(N_neg)  #in eV
            flag_negat=logical_or(en_eV_le<0., en_eV_le>4.)
            N_neg=sum(flag_negat);
        sigma_le=En_impact_eV[flag_low_energy]/4.
        en_eV_le=(en_eV_le+2.)*sigma_le
    
        
        #generate high energy
        en_eV_he=lognormal(mufit,sigmafit,N_high_ene);
        flag_above_th=logical_or(en_eV_he>E_th,(en_eV_he-En_impact_eV[flag_high_energy])>0)
        Nabove_th=sum(flag_above_th);
        
        while Nabove_th>0:
            en_eV_he[flag_above_th]=lognormal(mufit,sigmafit,Nabove_th);
            flag_above_th=logical_or(en_eV_he>E_th,(en_eV_he-En_impact_eV[flag_high_energy])>0)
            Nabove_th=sum(flag_above_th);
            
        en_eV[flag_high_energy]=en_eV_he
        en_eV[flag_low_energy]=en_eV_le
    
    return en_eV




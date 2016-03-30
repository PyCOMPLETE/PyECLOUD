#----------------------------------------------------------------------
#                                                                      
#                           CERN                                       
#                                                                      
#     European Organization for Nuclear Research                       
#                                                                      
#     
#     This file is part of the code:
#                                                                      		            
#		           PyECLOUD Version 5.0.2                     
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
#
#----------------------------------------------------------------------
#
#	SEY Model made by Jan Sopousek
#                                                
#----------------------------------------------------------------------

from numpy import sqrt, exp
from numpy.random import rand
import numpy as np
from operator import itemgetter


class SEY_model_FromFile:
    def __init__(self, FileName,pyecl_input_folder):
            self.FileName=FileName
            #The File Format:Energy[eV] Reflected Secondary
            #The Numbers Reflected and Secondary correspond to measured data by equation:
            #Measured=Reflected+(1-Reflected)*Secondary
            f=open(pyecl_input_folder+'/'+FileName)
            Data=[] 
            for line in f.readlines():
            	ls=line.split()
            	Energy=float(ls[0])
            	Reflected=float(ls[1])
            	Delta=float(ls[2])
            	
            	if Reflected == 1.:
					Secondary = 0.
            	else:
					Secondary = (Delta-Reflected)/(1-Reflected)
            	
            	Data.append([Energy,Reflected,Secondary])
            	
            Data=sorted(Data, key=itemgetter(0))
            self.CreateInterpolationMap(np.array(Data))
            
            
            print 'Secondary emission model: FromFile:',FileName,' No. points:',self.N
			
			
			
			
    def CreateInterpolationMap(self,Data):
        DeltaEs=Data[1:,0]-Data[0:-1,0]
        DeltaE=min(DeltaEs[DeltaEs>0])
		
        self.N=int(min(Data[Data.shape[0]-1,0]/DeltaE,100000)+1)	
        self.DeltaE=1.*Data[Data.shape[0]-1,0]/(self.N-1)
		
        Last=np.array([1,0])
        self.InterpolationMap=np.zeros((self.N,2))
        j=0
        i=0
        while i<self.N:
            while j<Data.shape[0]-1 and Data[j,0]<=self.DeltaE*i:
                j+=1			
            
            if j==0:
                
                self.InterpolationMap[i,:]=Data[0,1:]
            else:                   
                self.InterpolationMap[i,:]=(Data[j-1,1:]*(Data[j,0]-self.DeltaE*i)+Data[j,1:]*(self.DeltaE*i-Data[j-1,0]))/(Data[j,0]-Data[j-1,0])
            i+=1
				
    def GetR(self,E):
        i=np.floor(E/self.DeltaE).astype('int')
        MaskOver=(i>self.N-2)
        i=np.minimum(i,self.N-2)
        Results=(self.InterpolationMap[i,0]*(self.DeltaE*(i+1)-E)+self.InterpolationMap[i+1,0]*(E-self.DeltaE*(i)))/self.DeltaE
        Results[MaskOver]=((self.InterpolationMap[self.N-1,0])+0*Results)[MaskOver]
        return Results

    def GetS(self,E):
        i=np.floor(E/self.DeltaE).astype('int')
        MaskOver=(i>self.N-2)
        i=np.minimum(i,self.N-2)
        Results=(self.InterpolationMap[i,1]*(self.DeltaE*(i+1)-E)+self.InterpolationMap[i+1,1]*(E-self.DeltaE*(i)))/self.DeltaE
        Results[MaskOver]=((self.InterpolationMap[self.N-1,0])+0*Results)[MaskOver]
        return Results
	

    def Yield(self,E_impact_eV,costheta):
		costheta=np.abs(costheta)
		
		E=E_impact_eV/(1.+0.7*(1.-costheta))
		

		Secondary=self.GetS(E)*exp(0.5*(1.-costheta))

		Reflected=self.GetR(E_impact_eV)
		R=Reflected[Reflected<0]
		
		MaskTotal=(Reflected>=1)
		One=np.ones(len(E))
		n=One*2. #if set to 1 there is problem when costheta=0
		n[~MaskTotal]=(1.+sqrt(Reflected[~MaskTotal]))/(1.-sqrt(Reflected[~MaskTotal]))
		costhetat=sqrt(1.-(1.-costheta**2)/n**2)
		R=0.5*(((costheta-n*costhetat)/(costheta+n*costhetat))**2+((costhetat-n*costheta)/(costhetat+n*costheta))**2)
		#print costheta[R>1],costhetat[R>1],n[R>1]
		Reflected[~MaskTotal]=R[~MaskTotal]
		
		Reflected[Reflected>1]=Reflected[Reflected>1]*0.+1.
		
		delta=Secondary*(One-Reflected)+Reflected
		
		ref_frac=0.*delta
		mask_non_zero=(delta>0)
		ref_frac[mask_non_zero]=Reflected[mask_non_zero]/delta[mask_non_zero];

		return delta, ref_frac

 
    def SEY_process(self,nel_impact,E_impact_eV, costheta_impact, i_impact):
            yiel, ref_frac=self.Yield(E_impact_eV,costheta_impact);
            flag_elast=(rand(len(ref_frac))<ref_frac);
            flag_truesec=~(flag_elast);
            nel_emit=nel_impact*yiel;
            
            return  nel_emit, flag_elast, flag_truesec   

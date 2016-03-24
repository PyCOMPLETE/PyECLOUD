#----------------------------------------------------------------------
#                                                                      
#                           CERN                                       
#                                                                      
#     European Organization for Nuclear Research                       
#                                                                      
#     
#     This file is part of the code:
#                                                                      		            
#		           PyECLOUD Version 4.34                     
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
#	Rect-ellipse chamber made by Jan Sopousek
#                                                
#----------------------------------------------------------------------


from numpy import sum, sqrt, arctan2, sin, cos, isnan
import numpy as np

class rectellip_cham_geom_object:
	def __init__(self, filename_chm):
		self.a=filename_chm[0]
		self.b=filename_chm[1]
		self.c=filename_chm[2]
		self.d=filename_chm[3]
		
		self.x_aper=filename_chm[2]
		self.y_aper=filename_chm[3]
        
		self.N_mp_impact=0
		self.N_mp_corrected=0
		self.chamb_type='rectellip'


    
	def is_outside(self, x_mp, y_mp):
		A=np.logical_or(np.logical_or(np.logical_or(np.logical_or((((x_mp/self.c)**2 + (y_mp/self.d)**2)>=1),(x_mp<=-self.a)),(x_mp>=self.a)),(y_mp<=-self.b)),(y_mp>=self.b))
		return A



    
	def impact_point_and_normal(self, x_in, y_in, z_in, x_out, y_out, z_out, resc_fac=0.99, flag_robust=True):
	

    
		self.N_mp_impact=self.N_mp_impact+len(x_in)


		r=np.array([x_in,y_in])
		v=np.array([x_out-x_in,y_out-y_in])

		M=np.array([[1./self.c**2,0],[0,1./self.d**2]])

		A=np.dot(np.ones((1,2)),v*np.dot(M,v))
		B=np.dot(np.ones((1,2)),v*np.dot(M,r))+np.dot(np.ones((1,2)),r*np.dot(M,v))
		C=np.dot(np.ones((1,2)),r*np.dot(M,r))-np.ones((1,len(x_in)))



		t0=(-B+np.sqrt(B**2-4*A*C))/(2*A)
		t0=t0[0,:]	

			
		x_int=t0*x_out+(1-t0)*x_in;
		y_int=t0*y_out+(1-t0)*y_in;		
		z_int=t0*z_out+(1-t0)*z_in;	
			

		ZeroArray=0.*x_int
		Dx=0*ZeroArray
		Dy=0*ZeroArray
		
		MaskYZero=(y_int==0)
		Dx[MaskYZero]=0.
		Dy[MaskYZero]=1.
		Dx[~MaskYZero]=1.
		Dy[~MaskYZero]=(-x_int[~MaskYZero]*self.d**2/self.c**2/y_int[~MaskYZero])
		tt=0*t0
		
		MaskLeft=(x_in<self.a) & (self.a<x_out)
		tt[MaskLeft]=(self.a-x_in[MaskLeft])/(x_out[MaskLeft]-x_in[MaskLeft])
		MaskLeft[MaskLeft]=(tt[MaskLeft]<t0[MaskLeft])
		Dx[MaskLeft]=0.
		Dy[MaskLeft]=1.
		t0[MaskLeft]=tt[MaskLeft]

		
		MaskRight=(-self.a<x_in) & (x_out<-self.a)
		tt[MaskRight]=(-self.a-x_in[MaskRight])/(x_out[MaskRight]-x_in[MaskRight])
		MaskRight[MaskRight]=(tt[MaskRight]<t0[MaskRight])
		Dx[MaskRight]=0.
		Dy[MaskRight]=1.
		t0[MaskRight]=tt[MaskRight]

		
		MaskDown=(y_in<self.b) & (self.b<y_out)
		tt[MaskDown]=(self.b-y_in[MaskDown])/(y_out[MaskDown]-y_in[MaskDown])
		MaskDown[MaskDown]=(tt[MaskDown]<t0[MaskDown])
		Dx[MaskDown]=1.
		Dy[MaskDown]=0.
		t0[MaskDown]=tt[MaskDown]

		MaskUp=(-self.b<y_in) & (y_out<-self.b)
		tt[MaskUp]=(-self.b-y_in[MaskUp])/(y_out[MaskUp]-y_in[MaskUp])
		MaskUp[MaskUp]=(tt[MaskUp]<t0[MaskUp])
		Dx[MaskUp]=1.
		Dy[MaskUp]=0.
		t0[MaskUp]=tt[MaskUp]

			
		t0=resc_fac*t0;           
		t0[t0<1.e-2]=0;
			
		x_int=t0*x_out+(1-t0)*x_in;
		y_int=t0*y_out+(1-t0)*y_in;
		z_int=t0*z_out+(1-t0)*z_in;	




		Nx=-Dy;
		Ny=Dx;


		neg_flag=((Nx*x_int+Ny*y_int)>0);

		Nx[neg_flag]=-Nx[neg_flag];
		Ny[neg_flag]=-Ny[neg_flag];


		nor=sqrt(Nx*Nx+Ny*Ny);
		

		Nx=Nx/nor;
		Ny=Ny/nor;

		i_found=None

		return  x_int,y_int,z_int,Nx,Ny, i_found


#----------------------------------------------------------------------
#                                                                      
#                           CERN                                       
#                                                                      
#     European Organization for Nuclear Research                       
#                                                                      
#     
#     This file is part of the code:
#                                                                      		            
#		           PyECLOUD Version 4.23               
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



from numpy import sum, sqrt, arctan2, sin, cos, isnan

class ellip_cham_geom_object:
    def __init__(self, x_aper, y_aper, flag_verbose_file=True):
        self.x_aper = x_aper
        self.y_aper = y_aper
        
        self.N_mp_impact=0
        self.N_mp_corrected=0
        self.chamb_type='ellip'
        
        self.flag_verbose_file = flag_verbose_file

    
    def is_outside(self, x_mp, y_mp):
        return (((x_mp/self.x_aper)**2 + (y_mp/self.y_aper)**2)>=1);
    
    def impact_point_and_normal(self, x_in, y_in, z_in, x_out, y_out, z_out, resc_fac=0.99, flag_robust=True):
    
        self.N_mp_impact=self.N_mp_impact+len(x_in)
        
        a=self.x_aper
        b=self.y_aper
        
        x_insq=x_in*x_in;
        y_insq=y_in*y_in;
        x_outsq=x_out*x_out;
        y_outsq=y_out*y_out;
        y_in__y_out=y_in*y_out;
        x_in__x_out=x_in*x_out;
        
        
        t0=(a**2*y_insq - a**2*y_in__y_out + \
          sqrt(a**4*b**2*y_insq - 2*a**4*b**2*y_in__y_out \
              + a**4*b**2*y_outsq +  \
            a**2*b**4*x_insq - 2*a**2*b**4*x_in__x_out +  \
            a**2*b**4*x_outsq - a**2*b**2*x_insq* \
            y_outsq + 2*a**2*b**2*x_in__x_out*y_in__y_out -  \
            a**2*b**2*x_outsq*y_insq) + b**2*x_insq -  \
            b**2*x_in__x_out)/(a**2*y_insq -  \
            2*a**2*y_in__y_out + a**2*y_outsq + b**2*x_insq -  \
            2*b**2*x_in__x_out + b**2*x_outsq);
        
        
    
        # Handle pathological cases
        mask_nan = isnan(t0)
        if sum(mask_nan)>0:
            t0[mask_nan]=0.
            if self.flag_verbose_file:
                x_nan=x_in[mask_nan]
                y_nan=y_in[mask_nan]
                fbckt=open('bcktr_errors.txt','a')
                for ii_bk in xrange(len(y_nan)):
                    fbckt.write('%e\t%e\tnan\n'%(x_nan[ii_bk],y_nan[ii_bk]))
                fbckt.close()
        
        t0=resc_fac*t0;
           
        t0[t0<1.e-2]=0;
        
        flag_ident=(((x_in-x_out)*(x_in-x_out)+(y_in-y_out)*(y_in-y_out))/(x_out*x_out+y_out*y_out))<1e-8;
        t0[flag_ident]=0;
        
       
        
        if sum(abs(t0.imag))>0:
            print 'imag detected'
            raise ValueError('Backtracking: complex t0!!!!')
        
        x_int=t0*x_out+(1-t0)*x_in;
        y_int=t0*y_out+(1-t0)*y_in;
        z_int=t0*z_out+(1-t0)*z_in;
        
        
        if flag_robust:
			flag_impact=(((x_int/a)**2 + (y_int/b)**2)>1);
			
			
			if flag_impact.any():
				self.N_mp_corrected = self.N_mp_corrected + sum(flag_impact)
				ntrials=10
				while (sum(flag_impact)>0 and ntrials>0):
					t0[flag_impact]=0.9*t0[flag_impact];
					x_int=t0*x_out+(1-t0)*x_in;
					y_int=t0*y_out+(1-t0)*y_in;
					z_int=t0*z_out+(1-t0)*z_in;
					
					flag_impact=(((x_int/a)**2 + (y_int/b)**2)>1);
					ntrials=ntrials-1
				
			flag_impact=(((x_int/a)**2 + (y_int/b)**2)>=1);
			if sum(flag_impact)>0:
				
				x_int_pat = x_int[flag_impact]
				y_int_pat = y_int[flag_impact]
				
				if self.flag_verbose_file:
					fbckt=open('bcktr_errors.txt','a')
					for ii_bk in xrange(len(x_int_pat)):
						fbckt.write('%e\t%e\n'%(x_int_pat[ii_bk],y_int_pat[ii_bk]))
				
					fbckt.close()
			
				
				x_pr=x_int_pat/a
				y_pr=y_int_pat/b
				
				r_pr=sqrt(x_pr**2+y_pr**2)
				
				x_pr=0.99*x_pr/r_pr
				y_pr=0.99*y_pr/r_pr
				
				x_int[flag_impact] = x_pr*a
				y_int[flag_impact] = y_pr*b
				
				
				
				
			flag_impact=(((x_int/a)**2 + (y_int/b)**2)>=1);
			if sum(flag_impact)>0:   
				print 'err inside'
				raise ValueError('Outside after backtracking!!!!')
			
			
    
        
        
        par_cross=arctan2(a*y_int,b*x_int);
        
        Dx=-a*sin(par_cross);
        Dy=b*cos(par_cross);
        
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


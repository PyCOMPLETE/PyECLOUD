#----------------------------------------------------------------------
#                                                                      
#                           CERN                                       
#                                                                      
#     European Organization for Nuclear Research                       
#                                                                      
#     
#     This file is part of the code:
#                                                                      		            
#		           PyECLOUD Version 4.13                     
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

import numpy as np
import scipy.sparse as scsp
from scipy.sparse.linalg import spsolve
import rhocompute as rhocom
import int_field_for as iff
import scipy.sparse.linalg as ssl

class space_charge:
    
    def __init__(self,chamb, Dh, Dt_sc=None):
        
        print 'Start LU space charge init.'
      
        xg=np.arange(0,chamb.x_aper+2.*Dh,Dh,float)  
        xgr=xg[1:]
        xgr=xgr[::-1]#reverse array
        xg=np.concatenate((-xgr,xg),0)
        Nxg=len(xg);
        bias_x=min(xg);
        
        yg=np.arange(0,chamb.y_aper+2.*Dh,Dh,float)  
        ygr=yg[1:]
        ygr=ygr[::-1]#reverse array
        yg=np.concatenate((-ygr,yg),0)
        Nyg=len(yg);
        bias_y=min(yg);
        
        [xn, yn]=np.meshgrid(xg,yg)
        
        xn=xn.T
        xn=xn.flatten()
        
        yn=yn.T
        yn=yn.flatten()
        #% xn and yn are stored such that the external index is on x 
        
        flag_outside_n=chamb.is_outside(xn,yn)
        flag_inside_n=~(flag_outside_n)
        #flag_inside_n=(((xn/x_aper)**2 + (yn/y_aper)**2)<1);
        #flag_outside_n= ~(flag_inside_n);
        
        flag_outside_n_mat=np.reshape(flag_outside_n,(Nyg,Nxg),'F');
        flag_outside_n_mat=flag_outside_n_mat.T
        [gx,gy]=np.gradient(np.double(flag_outside_n_mat));
        gradmod=abs(gx)+abs(gy);
        flag_border_mat=np.logical_and((gradmod>0), flag_outside_n_mat);
        
        
        A=scsp.lil_matrix((Nxg*Nyg,Nxg*Nyg)); #allocate a sparse matrix
        
        
        #~ A.setdiag(ones((Nxg*Nyg,1)))
        #set regular stencil on internal nodes
        for u in range(0,Nxg*Nyg):
            if np.mod(u, Nxg*Nyg/20)==0:
                print ('Mat. assembly %.0f'%(float(u)/ float(Nxg*Nyg)*100)+"""%""")
            if flag_inside_n[u]:
                A[u,u] = -4/(Dh*Dh);    #phi(i,j)
                A[u,u-1]=1/(Dh*Dh);     #phi(i-1,j)nx
                A[u,u+1]=1/(Dh*Dh);     #phi(i+1,j)
                A[u,u-Nyg]=1/(Dh*Dh);    #phi(i,j-1)
                A[u,u+Nyg]=1/(Dh*Dh);    #phi(i,j+1)
            else:
                A[u,u]=1.

    
        #A=A.tocsr() #convert to csr format
        luobj = ssl.splu(A.tocsc())
        
    
        self.Dh=Dh
        self.xg = xg
        self.Nxg = Nxg
        self.bias_x = bias_x
        self.yg = yg
        self.Nyg = Nyg
        self.bias_y = bias_y
        self.xn = xn
        self.yn = yn
        self.flag_inside_n = flag_inside_n
        self.flag_outside_n = flag_outside_n
        self.flag_outside_n_mat = flag_outside_n_mat
        self.flag_border_mat = flag_border_mat
        self.A = A
        self.luobj = luobj
        
        
        self.rho = np.zeros((self.Nxg,self.Nyg));
        self.phi = np.zeros((self.Nxg,self.Nyg));
        self.efx = np.zeros((self.Nxg,self.Nyg));
        self.efy = np.zeros((self.Nxg,self.Nyg));
        
        self.Dt_sc = Dt_sc
        self.t_last_recom=0.;
        self.U_sc_eV_stp=0.;
        
        self.flag_decimate=(self.Dt_sc is not None)
        self.flag_recomputed_sc=False
        
        print 'Done space charge init.'
                        
    
    def recompute_spchg_efield(self, MP_e, t_curr=None):
        
        flag_recompute=True              
        if self.flag_decimate:
            flag_recompute = (t_curr - self.t_last_recom)>=self.Dt_sc
        
        if flag_recompute:
            self.t_last_recom = t_curr
            
            #print t_curr
            
            if (MP_e.N_mp>0):
            
                qe=1.602176565e-19;
                eps0=8.8541878176e-12;
                
                rho=rhocom.compute_sc_rho(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp],MP_e.nel_mp[0:MP_e.N_mp],
                                          self.bias_x,self.bias_y,self.Dh, self.Nxg, self.Nyg)
                #rho=compute_sc_rho(N_mp,x_mp, y_mp,nel_mp, bias_x, bias_y, Dh, Nxg, Nyg)
                    
                rho=-qe*rho/(self.Dh*self.Dh);
                
                b=-rho.flatten()/eps0;
                b[(~self.flag_inside_n)]=0; #boundary condition
                
                #phi=spsolve(self.A, b);
                phi = self.luobj.solve(b)
                
                
                U_sc_eV_stp=-0.5*eps0*sum(b*phi)*self.Dh*self.Dh/qe
                
                phi=np.reshape(phi,(self.Nxg,self.Nyg))
                #
                #
                efx=np.zeros((self.Nxg,self.Nyg));
                efy=np.zeros((self.Nxg,self.Nyg));
                efx[1:self.Nxg-1,:] = phi[0:self.Nxg-2,:] - phi[2:self.Nxg,:];  #central difference on internal nodes
                efy[:,1:self.Nyg-1] = phi[:,0:self.Nyg-2] - phi[:,2:self.Nyg];  #central difference on internal nodes
                #% phi_mat=reshape(phi,Nyg,Nxg).';
                #% flag_border=not(insidemat)&((abs(efx)+abs(efy))>0);
                efx[self.flag_border_mat]=efx[self.flag_border_mat]*2;
                efy[self.flag_border_mat]=efy[self.flag_border_mat]*2;
                #
                efx = efx / (2*self.Dh);    #divide grid size
                efy = efy / (2*self.Dh);
                
                   
                self.rho = rho
                self.b = b
                self.phi = phi
                self.efx = efx
                self.efy = efy
                self.U_sc_eV_stp = U_sc_eV_stp
                
                self.Dt_sc
        self.flag_recomputed_sc=flag_recompute
        
    def compute_spchg_efield_from_rho(self, rho, flag_verbose = True):
        
        qe=1.602176565e-19;
        eps0=8.8541878176e-12;
                
        b=-rho.flatten()/eps0;
        b[(~self.flag_inside_n)]=0; #boundary condition
        
        if flag_verbose:
            print 'Start Linear System Solution.'
            
        phi=spsolve(self.A, b);
        
        U_sc_eV_stp=-0.5*eps0*sum(b*phi)*self.Dh*self.Dh/qe
        
        phi=np.reshape(phi,(self.Nxg,self.Nyg))
        #
        #
        if flag_verbose:
            print 'Start field computation.'
        
        efx=np.zeros((self.Nxg,self.Nyg));
        efy=np.zeros((self.Nxg,self.Nyg));
        efx[1:self.Nxg-1,:] = phi[0:self.Nxg-2,:] - phi[2:self.Nxg,:];  #central difference on internal nodes
        efy[:,1:self.Nyg-1] = phi[:,0:self.Nyg-2] - phi[:,2:self.Nyg];  #central difference on internal nodes
        #% phi_mat=reshape(phi,Nyg,Nxg).';
        #% flag_border=not(insidemat)&((abs(efx)+abs(efy))>0);
        efx[self.flag_border_mat]=efx[self.flag_border_mat]*2;
        efy[self.flag_border_mat]=efy[self.flag_border_mat]*2;
        #
        efx = efx / (2*self.Dh);    #divide grid size
        efy = efy / (2*self.Dh);
        
           
        self.rho = rho
        self.b = b
        self.phi = phi
        self.efx = efx
        self.efy = efy
        self.U_sc_eV_stp = U_sc_eV_stp
                                
    def get_sc_eletric_field(self, MP_e):
            
        if MP_e.N_mp>0:    
            ## compute beam electric field
            Ex_sc_n, Ey_sc_n = iff.int_field(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp],self.bias_x,self.bias_y,self.Dh,
                                         self.Dh, self.efx, self.efy)
                       
        else:
            Ex_sc_n=0.
            Ey_sc_n=0.
            
        return Ex_sc_n, Ey_sc_n
        
          

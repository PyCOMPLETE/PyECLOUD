#----------------------------------------------------------------------
#                                                                      
#                           CERN                                       
#                                                                      
#     European Organization for Nuclear Research                       
#                                                                      
#     
#     This file is part of the code:
#                                                                      		    
# 
#		           PyECLOUD Version 4.27                     
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
from vectsum import vectsum

na = lambda x:np.array([x])

class space_charge:
    #@profile
    def __init__(self,chamb, Dh, Dt_sc=None, sparse_solver = 'scipy_slu'):
        
		print 'Start LU space charge init.'
		print 'Using Shortley-Weller boundary approx.'
		xg=np.arange(0,chamb.x_aper+5.*Dh,Dh,float)  
		xgr=xg[1:]
		xgr=xgr[::-1]#reverse array
		xg=np.concatenate((-xgr,xg),0)
		Nxg=len(xg);
		bias_x=min(xg);

		yg=np.arange(0,chamb.y_aper+4.*Dh,Dh,float)  
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

		flag_outside_n_mat=np.reshape(flag_outside_n,(Nyg,Nxg),'F');
		flag_outside_n_mat=flag_outside_n_mat.T
		[gx,gy]=np.gradient(np.double(flag_outside_n_mat));
		gradmod=abs(gx)+abs(gy);
		flag_border_mat=np.logical_and((gradmod>0), flag_outside_n_mat);
		flag_border_n = flag_border_mat.flatten()

		A=scsp.lil_matrix((Nxg*Nyg,Nxg*Nyg)); #allocate a sparse matrix
		Dx=scsp.lil_matrix((Nxg*Nyg,Nxg*Nyg)); #allocate a sparse matrix
		Dy=scsp.lil_matrix((Nxg*Nyg,Nxg*Nyg)); #allocate a sparse matrix

		list_internal_force_zero = []

		# Build A Dx Dy matrices 
		for u in range(0,Nxg*Nyg):
			if np.mod(u, Nxg*Nyg/20)==0:
				print ('Mat. assembly %.0f'%(float(u)/ float(Nxg*Nyg)*100)+"""%""")
			if flag_inside_n[u]:
				
				#Compute Shortley-Weller coefficients
				if flag_inside_n[u-1]: #phi(i-1,j)
					hw = Dh
				else:
					x_int,y_int,z_int,Nx_int,Ny_int, i_found_int = chamb.impact_point_and_normal(na(xn[u]), na(yn[u]), na(0.), na(xn[u-1]), na(yn[u-1]), na(0.), resc_fac=.995, flag_robust=False)
					hw = np.abs(y_int[0]-yn[u])
					
				if flag_inside_n[u+1]: #phi(i+1,j)
					he = Dh
				else:
					x_int,y_int,z_int,Nx_int,Ny_int, i_found_int = chamb.impact_point_and_normal(na(xn[u]), na(yn[u]), na(0.), na(xn[u+1]), na(yn[u+1]), na(0.), resc_fac=.995, flag_robust=False)
					he = np.abs(y_int[0]-yn[u])
				
				if flag_inside_n[u-Nyg]: #phi(i,j-1)
					hs = Dh
				else:
					x_int,y_int,z_int,Nx_int,Ny_int, i_found_int = chamb.impact_point_and_normal(na(xn[u]), na(yn[u]), na(0.), na(xn[u-Nyg]), na(yn[u-Nyg]), na(0.), resc_fac=.995, flag_robust=False)
					hs = np.abs(x_int[0]-xn[u])
					#~ print hs
				
				if flag_inside_n[u+Nyg]: #phi(i,j+1)
					hn = Dh
				else:
					x_int,y_int,z_int,Nx_int,Ny_int, i_found_int = chamb.impact_point_and_normal(na(xn[u]), na(yn[u]), na(0.), na(xn[u+Nyg]), na(yn[u+Nyg]), na(0.), resc_fac=.995, flag_robust=False)
					hn = np.abs(x_int[0]-xn[u])
					#~ print hn
				
				
				# Build A matrix
				if hn<Dh/100. or hs<Dh/100. or hw<Dh/100. or he<Dh/100.: # nodes very close to the bounday
					A[u,u] =1.
					list_internal_force_zero.append(u)
					#print u, xn[u], yn[u]
				else:
					A[u,u] = -(2./(he*hw)+2/(hs*hn))
					A[u,u-1]=2./(hw*(hw+he));     #phi(i-1,j)nx
					A[u,u+1]=2./(he*(hw+he));     #phi(i+1,j)
					A[u,u-Nyg]=2./(hs*(hs+hn));    #phi(i,j-1)
					A[u,u+Nyg]=2./(hn*(hs+hn));    #phi(i,j+1)
				
				
				
				# Build Dx matrix
				if hn<Dh/100.:
					if hs>=Dh/100.:
						Dx[u,u] = -1./hs
						Dx[u,u-Nyg]=1./hs
				elif hs<Dh/100.:
					if hn>=Dh/100.:
						Dx[u,u] = 1./hn
						Dx[u,u+Nyg]=-1./hn
				else:
					Dx[u,u] = (1./(2*hn)-1./(2*hs))
					Dx[u,u-Nyg]=1./(2*hs)
					Dx[u,u+Nyg]=-1./(2*hn)
					
					
				# Build Dy matrix	
				if he<Dh/100.:
					if hw>=Dh/100.:
						Dy[u,u] = -1./hw
						Dy[u,u-1]=1./hw
				elif hw<Dh/100.:
					if he>=Dh/100.:
						Dy[u,u] = 1./he
						Dy[u,u+1]=-1./(he)
				else:
					Dy[u,u] = (1./(2*he)-1./(2*hw))
					Dy[u,u-1]=1./(2*hw)
					Dy[u,u+1]=-1./(2*he)				

					
			else:
				# external nodes
				A[u,u]=1.
				if flag_border_n[u]:
					handle_border(u, flag_inside_n, Nxg, Nyg, xn, yn, chamb, Dh, Dx, Dy)	

		for u in list_internal_force_zero:
			handle_border(u, flag_inside_n, Nxg, Nyg, xn, yn, chamb, Dh, Dx, Dy)

		#~ A = A.tocsc()
		#~ Dx = Dx.tocsc()
		#~ Dy = Dy.tocsc()

		flag_force_zero = flag_outside_n.copy()	
		for ind in 	list_internal_force_zero:
			flag_force_zero[ind] = True
			
		flag_force_zero_mat=np.reshape(flag_force_zero,(Nyg,Nxg),'F');
		flag_force_zero_mat=flag_force_zero_mat.T
		[gxc,gyc]=np.gradient(np.double(flag_force_zero_mat));
		gradmodc=abs(gxc)+abs(gyc);
		flag_border_mat_c=np.logical_and((gradmodc>0), flag_force_zero_mat);
			
		sumcurr = np.sum(flag_border_mat_c, axis=0)
		jj_max_border = np.max((np.where(sumcurr>0))[0])
		jj_min_border = np.min((np.where(sumcurr>0))[0])

		sumcurr = np.sum(flag_border_mat_c, axis=1)# corrected in version 4.05. I it was: sumcurr = np.sum(flag_border_mat_c, axis=1)
		ii_max_border = np.max((np.where(sumcurr>0))[0])
		ii_min_border = np.min((np.where(sumcurr>0))[0])
			
		print 'Internal nodes with 0 potential'
		print list_internal_force_zero

		A=A.tocsr() #convert to csr format
		
		#Remove trivial equtions 
		diagonal = A.diagonal()
		N_full = len(diagonal)
		indices_non_id = np.where(diagonal!=1.)[0]
		N_sel = len(indices_non_id)

		Msel = scsp.lil_matrix((N_full, N_sel))
		for ii, ind in enumerate(indices_non_id):
			Msel[ind, ii] =1.
			
		Msel = Msel.tocsc()

		Asel = Msel.T*A*Msel
		Asel=Asel.tocsc()
		

		if sparse_solver == 'scipy_slu':
			print "Using scipy superlu solver..."
			luobj = ssl.splu(Asel.tocsc())
		elif sparse_solver == 'klu':
			print "Using klu solver..."
			try:
				import klulib.klu as klu
				luobj = klu.Klu(Asel.tocsc())
			except StandardError, e: 
				print "Got exception: ", e
				print "Falling back on scipy superlu solver:"
				luobj = ssl.splu(Asel.tocsc())
		else:
			raise ValueError('Solver not recognized!!!!\nsparse_solver must be "scipy_slu" or "klu"\n')
			
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
		self.flag_inside_n_mat = np.logical_not(flag_outside_n_mat)
		self.flag_border_mat = flag_border_mat
		self.flag_force_zero = flag_force_zero
		self.Asel = Asel
		self.luobj = luobj
		self.Dx = Dx.tocsc()
		self.Dy = Dy.tocsc()

		self.ii_max_border = ii_max_border
		self.ii_min_border = ii_min_border
		self.jj_max_border = jj_max_border
		self.jj_min_border = jj_min_border

		self.rho = np.zeros((self.Nxg,self.Nyg));
		self.phi = np.zeros((self.Nxg,self.Nyg));
		self.efx = np.zeros((self.Nxg,self.Nyg));
		self.efy = np.zeros((self.Nxg,self.Nyg));

		self.Dt_sc = Dt_sc
		self.t_last_recom=0.;
		self.U_sc_eV_stp=0.;

		self.flag_decimate=(self.Dt_sc is not None)
		self.flag_recomputed_sc=False
		
		self.Msel = Msel.tocsc()
		self.Msel_T = (Msel.T).tocsc()

		print 'Done space charge init.'
                        
    #@profile
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

                rho=-qe*rho/(self.Dh*self.Dh);
                
                self.compute_spchg_efield_from_rho(rho, flag_verbose = False)
                
                
        self.flag_recomputed_sc=flag_recompute

    #@profile    
    def compute_spchg_efield_from_rho(self, rho, flag_verbose = True):
        
        qe=1.602176565e-19;
        eps0=8.8541878176e-12;
                
        
        b=-rho.flatten()/eps0;
        b[(self.flag_force_zero)]=0; #boundary condition
		
        if flag_verbose:
            print 'Start Linear System Solution.'
        b_sel = self.Msel_T*b
        phi_sel = self.luobj.solve(b_sel)
        phi = self.Msel*phi_sel
		
        U_sc_eV_stp = -0.5*eps0*np.sum(b*phi)*self.Dh*self.Dh/qe
		
        if flag_verbose:
            print 'Start field computation.'
		
        efx = self.Dx*phi
        efy = self.Dy*phi
        phi=np.reshape(phi,(self.Nxg,self.Nyg))
        efx=np.reshape(efx,(self.Nxg,self.Nyg))
        efy=np.reshape(efy,(self.Nxg,self.Nyg))
        
        for jj in xrange(self.jj_max_border, self.Nyg):
			efx[:, jj]=efx[:, self.jj_max_border-1] 
			
        for jj in xrange(0, self.jj_min_border+1):
			efx[:, jj]=efx[:, self.jj_min_border+1] 
			
        for ii in xrange(self.ii_max_border, self.Nxg):
			efy[ii, :]=efy[self.ii_max_border-1, :] 
			
        for ii in xrange(0, self.ii_min_border+1):
			efy[ii,:]=efy[self.ii_min_border+1,:] 
		   
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

def handle_border(u, flag_inside_n, Nxg, Nyg, xn, yn, chamb, Dh, Dx, Dy):
	#print u
	jjj = np.floor(u/Nyg)
	
	if flag_inside_n[u+Nyg]: 
		if not flag_inside_n[u]:
			x_int,y_int,z_int,Nx_int,Ny_int, i_found_int = chamb.impact_point_and_normal( na(xn[u+Nyg]), na(yn[u+Nyg]), na(0.), 
				na(xn[u]), na(yn[u]), na(0.), resc_fac=.995, flag_robust=False)
			hs = np.abs(x_int[0]-xn[u+Nyg])
		else: #this is the case for internal nodes with zero potential (very close to the boundary)
			hs = Dh
		
		hn = Dh

		if hs<Dh/100.:
			Dx[u,u+Nyg] = (1./(hn))
			Dx[u,u+Nyg+Nyg]=-1./(hn)
			
			nnn=1
			while u-nnn*Nyg>=0:
				Dx[u-nnn*Nyg,u+Nyg] = (1./(hn))
				Dx[u-nnn*Nyg,u+Nyg+Nyg]=-1./(hn)
				nnn+=1
			
		else:
			Dx[u,u+Nyg] = (1./(2*hn)-1./(2*hs))
			Dx[u,u-Nyg+Nyg] = 1./(2*hs)
			Dx[u,u+Nyg+Nyg] = -1./(2*hn)
			
			nnn=1
			while u-nnn*Nyg>=0:
				Dx[u-nnn*Nyg,u+Nyg] = Dx[u,u+Nyg]
				Dx[u-nnn*Nyg,u-Nyg+Nyg] = Dx[u,u-Nyg+Nyg]
				Dx[u-nnn*Nyg,u+Nyg+Nyg] = Dx[u,u+Nyg+Nyg]
				nnn+=1

	
	elif flag_inside_n[u-Nyg]: 
		if not flag_inside_n[u]:
			x_int,y_int,z_int,Nx_int,Ny_int, i_found_int = chamb.impact_point_and_normal( na(xn[u-Nyg]), na(yn[u-Nyg]), na(0.), 
				na(xn[u]), na(yn[u]), na(0.), resc_fac=.995, flag_robust=False)
			hn = np.abs(x_int[0]-xn[u-Nyg])
		else:#this is the case for internal nodes with zero potential (very close to the boundary)
			hn = Dh
		
		hs = Dh

		if hn<Dh/100.:
			Dx[u,u-Nyg] = -1./(hs)
			Dx[u,u-Nyg-Nyg]=1./(hs)
			
			nnn=1
			while u+nnn*Nyg<Nxg*Nyg:
				Dx[u+nnn*Nyg,u-Nyg] = -1./(hs)
				Dx[u+nnn*Nyg,u-Nyg-Nyg]=1./(hs)
				nnn+=1
			
		else:
			Dx[u,u-Nyg] = (1./(2*hn)-1./(2*hs))
			Dx[u,u-Nyg-Nyg]=1./(2*hs)
			Dx[u,u+Nyg-Nyg]=-1./(2*hn)
			
			nnn=1
			while u+nnn*Nyg<Nxg*Nyg:
				Dx[u+nnn*Nyg,u-Nyg] = Dx[u,u-Nyg]
				Dx[u+nnn*Nyg,u-Nyg-Nyg] = Dx[u,u-Nyg-Nyg]
				Dx[u+nnn*Nyg,u+Nyg-Nyg] = Dx[u,u+Nyg-Nyg]
				nnn+=1
			
	if flag_inside_n[u+1]:
		if not flag_inside_n[u]: 
			x_int,y_int,z_int,Nx_int,Ny_int, i_found_int = chamb.impact_point_and_normal( na(xn[u+1]), na(yn[u+1]), na(0.), 
				na(xn[u]), na(yn[u]),na(0.), resc_fac=.995, flag_robust=False)
			hw = np.abs(y_int[0]-yn[u+1])
		else:#this is the case for internal nodes with zero potential (very close to the boundary)
			hw = Dh
			
		he = Dh
		
		if hw<Dh/100.:
			Dy[u,u+1] = (1./(he))
			Dy[u,u+1+1]=-1./(he)
			
			nnn=1
			while u-nnn>=(jjj)*Nyg:
				Dy[u-nnn*1,u+1] = (1./(he))
				Dy[u-nnn*1,u+1+1]=-1./(he)
				nnn+=1
		else:
			Dy[u,u+1] = (1./(2*he)-1./(2*hw))
			Dy[u,u-1+1] = 1./(2*hw)
			Dy[u,u+1+1] = -1./(2*he)
			
			nnn=1
			while u-nnn>=(jjj)*Nyg:
				#print nnn
				Dy[u-nnn,u+1] = Dy[u,u+1]
				Dy[u-nnn,u-1+1] = Dy[u,u-1+1]
				Dy[u-nnn,u+1+1] = Dy[u,u+1+1]
				nnn += 1
				
	elif flag_inside_n[u-1]: 
		if not flag_inside_n[u]:
			x_int,y_int,z_int,Nx_int,Ny_int, i_found_int = chamb.impact_point_and_normal( na(xn[u-1]), na(yn[u-1]), na(0.), 
				na(xn[u]), na(yn[u]),  na(0.), resc_fac=.995, flag_robust=False)
			he = np.abs(y_int[0]-yn[u-1])
		else:#this is the case for internal nodes with zero potential (very close to the boundary)
			he=Dh
			
		hw = Dh

		if he<Dh/100.:
			Dy[u,u-1] = -1./(hw)
			Dy[u,u-1-1]=1./(hw)
			
			nnn=1
			while u+nnn<(jjj+1)*Nyg:
				Dy[u+nnn,u-1] = -1./(hw)
				Dy[u+nnn,u-1-1]=1./(hw)
				nnn+=1
			
		else:
			Dy[u,u-1] = (1./(2*he)-1./(2*hw))
			Dy[u,u-1-1]=1./(2*hw)
			Dy[u,u+1-1]=-1./(2*he)
			
			nnn=1
			while u+nnn<(jjj+1)*Nyg:
				Dy[u+nnn,u-1] = Dy[u,u-1]
				Dy[u+nnn,u-1-1] = Dy[u,u-1-1]
				Dy[u+nnn,u+1-1] = Dy[u,u+1-1]
				nnn+=1
	return Dx, Dy

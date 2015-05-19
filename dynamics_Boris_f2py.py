#----------------------------------------------------------------------
#                                                                      
#                           CERN                                       
#                                                                      
#     European Organization for Nuclear Research                       
#                                                                      
#     
#     This file is part of the code:
#                                                                      		            
#		           PyECLOUD Version 4.25                     
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
import int_field_for as iff
from boris_step import boris_step



def crprod(bx, by, bz, cx,cy,cz):
    ax = by*cz-bz*cy
    ay = bz*cx-bx*cz
    az = bx*cy-by*cx
    
    return ax, ay, az


class B_none():
    
    def __init__(self, B0x, B0y, B0z):
        self.B0x=B0x
        self.B0y=B0y
        self.B0z=B0z
        
    def get_B(self, xn,yn):
            Bx_n = 0*xn + self.B0x
            By_n = 0*xn+ self.B0y
            Bz_n = 0*xn + self.B0z
            return Bx_n, By_n, Bz_n
        

class B_quad():
    
    def __init__(self, B0x, B0y, B0z, fact_Bmap):
        self.B0x=B0x
        self.B0y=B0y
        self.B0z=B0z
        self.fact_Bmap = fact_Bmap
        
        
    def get_B(self,xn,yn):
        Bx_n = self.fact_Bmap*yn.copy()
        By_n = self.fact_Bmap*xn.copy()
        Bx_n = Bx_n + self.B0x
        By_n = By_n + self.B0y
        Bz_n = 0*xn + self.B0z
        return Bx_n, By_n, Bz_n   
    
class B_file():
    
    def __init__(self, B0x, B0y, B0z, fact_Bmap, B_map_file):
        print 'Loading B map'
        dict_Bmap=sio.loadmat(B_map_file)

        self.Bmap_x = fact_Bmap*squeeze(dict_Bmap['Bx'].real)
        self.Bmap_y = fact_Bmap*squeeze(dict_Bmap['By'].real)
        self.xx=squeeze(dict_Bmap['xx'].T)
        self.yy=squeeze(dict_Bmap['yy'].T)
        self.B0x=B0x
        self.B0y=B0y
        self.B0z=B0z 
                 
        self.xmin=min(self.xx);
        self.ymin=min(self.yy);
        self.dx=self.xx[1]-self.xx[0];
        self.dy=self.yy[1]-self.yy[0];
    
    def get_B(self, xn,yn):
                Bx_n,By_n = iff.int_field(xn,yn,self.xmin,self.ymin,\
                                      self.dx,self.dy,self.Bmap_x,self.Bmap_y)
                # the rescaling factor has already been applied to the map
                Bx_n = Bx_n + self.B0x
                By_n = By_n + self.B0y
                Bz_n = 0*xn + self.B0z
                return Bx_n, By_n, Bz_n
            
    


class pusher_Boris():
    
    def __init__(self, Dt, B0x, B0y, B0z, \
                 B_map_file, fact_Bmap, Bz_map_file, N_sub_steps=1):
        
        print "Tracker: Boris"
        
        self.N_sub_steps = N_sub_steps
        self.Dtt = Dt / float(N_sub_steps)
        self.B0x = B0x
        self.B0y = B0y
        self.B0z = B0z
        self.fact_Bmap = fact_Bmap

        
        if B_map_file is None:
            self.B_ob = B_none(B0x, B0y, B0z)
            
        elif B_map_file is 'analytic_qaudrupole_unit_grad':
            print "B map analytic quadrupole"
            self.B_ob = B_quad(B0x, B0y, B0z, fact_Bmap)

        else:
            self.B_ob = B_file(B0x, B0y, B0z, fact_Bmap, B_map_file)
            
              
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
        
            if  Ez_n==0.:
                Ez_n = 0.*xn1
                
            for ii in range(self.N_sub_steps):
                Bx_n, By_n, Bz_n = self.B_ob.get_B(xn1,yn1)
                boris_step(self.Dtt,xn1,yn1,zn1,vxn1,vyn1,vzn1,
                           Ex_n,Ey_n,Ez_n,Bx_n,By_n,Bz_n)
                
            
        
            MP_e.x_mp[0:MP_e.N_mp] = xn1
            MP_e.y_mp[0:MP_e.N_mp] = yn1
            MP_e.z_mp[0:MP_e.N_mp] = zn1
            MP_e.vx_mp[0:MP_e.N_mp] = vxn1
            MP_e.vy_mp[0:MP_e.N_mp] = vyn1
            MP_e.vz_mp[0:MP_e.N_mp]  = vzn1      
                
            
            
        return MP_e
        
    def stepcustomDt(self, MP_e, Ex_n,Ey_n, Ez_n=0., Dt_substep=None, N_sub_steps=None):
        
        
        if MP_e.N_mp>0:
            
            xn1 = MP_e.x_mp[0:MP_e.N_mp]
            yn1 = MP_e.y_mp[0:MP_e.N_mp]
            zn1 = MP_e.z_mp[0:MP_e.N_mp]
            vxn1 = MP_e.vx_mp[0:MP_e.N_mp]
            vyn1 = MP_e.vy_mp[0:MP_e.N_mp]
            vzn1 = MP_e.vz_mp[0:MP_e.N_mp]
        
            if  Ez_n==0.:
                Ez_n = 0.*xn1
                
            for ii in range(N_sub_steps):
                Bx_n, By_n, Bz_n = self.B_ob.get_B(xn1,yn1)
                boris_step(Dt_substep,xn1,yn1,zn1,vxn1,vyn1,vzn1,
                           Ex_n,Ey_n,Ez_n,Bx_n,By_n,Bz_n)
                
            
        
            MP_e.x_mp[0:MP_e.N_mp] = xn1
            MP_e.y_mp[0:MP_e.N_mp] = yn1
            MP_e.z_mp[0:MP_e.N_mp] = zn1
            MP_e.vx_mp[0:MP_e.N_mp] = vxn1
            MP_e.vy_mp[0:MP_e.N_mp] = vyn1
            MP_e.vz_mp[0:MP_e.N_mp]  = vzn1      
                
            
            
        return MP_e       
     

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
#                   PyECLOUD Version 6.7.0
#
#
#     Main author:          Giovanni IADAROLA
#                           BE-ABP Group
#                           CERN
#                           CH-1211 GENEVA 23
#                           SWITZERLAND
#                           giovanni.iadarola@cern.ch
#
#     Author list:          Eleanora Belli
#                           Philipp Dijkstal
#                           Lotta Mether
#                           Annalisa Romano
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

na = lambda x:np.array([x])

class space_charge:
    #@profile

    def __init__(self,chamb, Dh, Dt_sc=None, PyPICmode = 'FiniteDifferences_ShortleyWeller' ,sparse_solver = 'scipy_slu',
                    f_telescope=None, target_grid=None, N_nodes_discard=None, N_min_Dh_main=None):

        print 'Start space charge init.'

        if PyPICmode == 'FiniteDifferences_ShortleyWeller':
            import PyPIC.FiniteDifferences_ShortleyWeller_SquareGrid as PIC_FDSW
            self.PyPICobj = PIC_FDSW.FiniteDifferences_ShortleyWeller_SquareGrid(chamb = chamb, Dh = Dh, sparse_solver = sparse_solver)
            #To be replaced by a property to make it general (from PyPIC modules not having xn, yn)
            self.xn = self.PyPICobj.xn
            self.yn = self.PyPICobj.yn
        elif PyPICmode == 'ShortleyWeller_WithTelescopicGrids':
            import PyPIC.FiniteDifferences_ShortleyWeller_SquareGrid as PIC_FDSW
            PyPICmain = PIC_FDSW.FiniteDifferences_ShortleyWeller_SquareGrid(chamb = chamb, Dh = Dh, sparse_solver = sparse_solver)
            import PyPIC.MultiGrid as PIC_MG
            self.PyPICobj = PIC_MG.AddTelescopicGrids(pic_main = PyPICmain, f_telescope = f_telescope, target_grid = target_grid,
                                        N_nodes_discard = N_nodes_discard, N_min_Dh_main = N_min_Dh_main, sparse_solver = sparse_solver)
            self.xn = None #not implemented in this mode (for now)
            self.yn = None #not implemented in this mode (for now)

        elif PyPICmode == 'FiniteDifferences_Staircase':
            import PyPIC.FiniteDifferences_Staircase_SquareGrid as PIC_FDSQ
            self.PyPICobj = PIC_FDSQ.FiniteDifferences_Staircase_SquareGrid(chamb = chamb, Dh = Dh, sparse_solver = sparse_solver)
            #To be replaced by a property to make it general (from PyPIC modules not having xn, yn)
            self.xn = self.PyPICobj.xn
            self.yn = self.PyPICobj.yn
        elif PyPICmode == 'FFT_PEC_Boundary':
            if chamb.chamb_type != 'rect':
                raise ValueError('''PyPICmode = 'FFT_PEC_Boundary' can be used only if chamb_type = 'rect' ''' )
            import PyPIC.FFT_PEC_Boundary_SquareGrid as PIC_FFT_PEC
            self.PyPICobj = PIC_FFT_PEC.FFT_PEC_Boundary_SquareGrid(x_aper = chamb.x_aper, y_aper = chamb.y_aper, Dh = Dh)
            #To be replaced by a property to make it general (from PyPIC modules not having xn, yn)
            self.xn = None #not implemented in this mode (for now)
            self.yn = None #not implemented in this mode (for now)
        elif PyPICmode == 'FFT_OpenBoundary':
            if chamb.chamb_type != 'rect':
                raise ValueError('''PyPICmode = 'FFT_OpenBoundary' can be used only if chamb_type = 'rect' ''' )
            import PyPIC.FFT_OpenBoundary as PIC_FFT_Open
            if len(np.atleast_1d(Dh)) == 2:
                self.PyPICobj = PIC_FFT_Open.FFT_OpenBoundary(x_aper = chamb.x_aper, y_aper = chamb.y_aper, dx = Dh[0], dy = Dh[1])
            else:
                self.PyPICobj = PIC_FFT_Open.FFT_OpenBoundary(x_aper = chamb.x_aper, y_aper = chamb.y_aper, Dh = Dh)

            #To be replaced by a property to make it general (from PyPIC modules not having xn, yn)
            self.xn = None #not implemented in this mode (for now)
            self.yn = None #not implemented in this mode (for now)
        else:
            raise ValueError('PyPICmode not racognized')

        self.Dh=self.PyPICobj.Dh
        self.xg = self.PyPICobj.xg
        self.Nxg = self.PyPICobj.Nxg
        self.bias_x = self.PyPICobj.bias_x
        self.yg = self.PyPICobj.yg
        self.Nyg = self.PyPICobj.Nyg
        self.bias_y = self.PyPICobj.bias_y


        self.Dt_sc = Dt_sc
        self.t_last_recom=0.;

        self.U_sc_eV_stp = 0.


        self.flag_decimate=(self.Dt_sc is not None)
        self.flag_recomputed_sc=False

        print 'Done space charge init.'

    @property
    def rho(self):
        return self.PyPICobj.rho

    @property
    def phi(self):
        return self.PyPICobj.phi

    @property
    def efx(self):
        return self.PyPICobj.efx

    @property
    def efy(self):
        return self.PyPICobj.efy


    #@profile
    def recompute_spchg_efield(self, MP_e, t_curr=None, force=False):

        flag_recompute=True
        if self.flag_decimate:
            flag_recompute = (t_curr - self.t_last_recom)>=self.Dt_sc

        if flag_recompute or force:
            self.t_last_recom = t_curr

            self.PyPICobj.scatter_and_solve(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp],MP_e.nel_mp[0:MP_e.N_mp], charge = MP_e.charge)
            #~ U_sc_eV_stp = -0.5*eps0*np.sum(b*phi)*self.Dh*self.Dh/qe

        self.flag_recomputed_sc=flag_recompute


    #@profile
    def compute_spchg_efield_from_rho(self, rho, flag_verbose = True):
        self.PyPICobj.solve(rho = rho, flag_verbose = flag_verbose)


    def get_sc_eletric_field(self, MP_e):
        Ex_sc_n, Ey_sc_n = self.PyPICobj.gather(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp])
        return Ex_sc_n, Ey_sc_n


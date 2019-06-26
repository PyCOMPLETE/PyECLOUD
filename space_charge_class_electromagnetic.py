#-Begin-preamble-------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 7.7.1
#
#
#     Main author:          Giovanni IADAROLA
#                           BE-ABP Group
#                           CERN
#                           CH-1211 GENEVA 23
#                           SWITZERLAND
#                           giovanni.iadarola@cern.ch
#
#     Contributors:         Eleonora Belli
#                           Philipp Dijkstal
#                           Lotta Mether
#                           Annalisa Romano
#                           Giovanni Rumolo
#                           Eric Wulff
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
#-End-preamble---------------------------------------------------------

import numpy as np
from space_charge_class import space_charge
from scipy.constants import c, epsilon_0, mu_0
import int_field_for as iff

na = lambda x: np.array([x])

class space_charge_electromagnetic(space_charge, object):

    def __init__(self, chamb, Dh, Dt_sc=None, PyPICmode='FiniteDifferences_ShortleyWeller' , sparse_solver='scipy_slu',
                 f_telescope=None, target_grid=None, N_nodes_discard=None, N_min_Dh_main=None, gamma = 479.):

        super(space_charge_electromagnetic, self).__init__(chamb, Dh, Dt_sc, PyPICmode , sparse_solver,
                     f_telescope, target_grid, N_nodes_discard, N_min_Dh_main)

        self.state_Ax = self.PyPICobj.get_state_object()
        self.state_Ay = self.PyPICobj.get_state_object()

        # Ax and Ay at previous steps for computation of dAx_dz and dAy_dz
        self.Ax_old_grid = np.zeros((self.Nxg,self.Nyg))
        self.Ay_old_grid = np.zeros((self.Nxg,self.Nyg))
        self.dAx_grid_dt = np.zeros((self.Nxg,self.Nyg))
        self.dAy_grid_dt = np.zeros((self.Nxg,self.Nyg))

        self.gamma = gamma
        self.beta = np.sqrt(1-1/(gamma*gamma))

    def recompute_spchg_efield(self, MP_e, flag_solve=True, flag_reset=True):
        # scatter rho
        self.PyPICobj.scatter(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp], MP_e.nel_mp[0:MP_e.N_mp],
                charge=MP_e.charge, flag_add=not(flag_reset))

        # scatter currents
        self.state_Ax.scatter(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp],
                epsilon_0 * mu_0 * MP_e.nel_mp[0:MP_e.N_mp] * MP_e.vx_mp[0:MP_e.N_mp],
                charge=MP_e.charge, flag_add=not(flag_reset))
        self.state_Ay.scatter(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp],
                epsilon_0 * mu_0 * MP_e.nel_mp[0:MP_e.N_mp] * MP_e.vy_mp[0:MP_e.N_mp],
                charge=MP_e.charge, flag_add=not(flag_reset))

        # solve
        if flag_solve:
            self.PyPICobj.solve()
            self.PyPICobj.solve_states([self.state_Ax, self.state_Ay])

        Ax_grid = self.state_Ax.phi
        Ay_grid = self.state_Ax.phi

        #compute time derivatives
        self.dAx_grid_dt = (Ax_grid -  self.Ax_old_grid)/self.Dt_sc
        self.dAy_grid_dt = (Ax_grid -  self.Ax_old_grid)/self.Dt_sc


        self.Ax_old_grid = Ax_grid
        self.Ay_old_grid = Ay_grid

    def get_sc_em_field(self, MP_e):
        #compute un-primed potentials
        _, dAx_dy = self.state_Ax.gather(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp])
        dAy_dx, _ = self.state_Ay.gather(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp])
        dphi_dx, dphi_dy = self.PyPICobj.gather(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp])

        #fix signs and make primed
        dphi_prime_dx = -self.gamma*dphi_dx
        dphi_prime_dy = -self.gamma*dphi_dy
        dAx_prime_dy = -self.gamma*dAx_dy
        dAy_prime_dx = -self.gamma*dAy_dx
        #compute longitudinal magnetic potential
        dAs_prime_dx = -self.beta/c*dphi_prime_dx
        dAs_prime_dy = -self.beta/c*dphi_prime_dy
        # interpolate time derivatives to the particle positions
        dAx_dt, dAy_dt = iff.int_field(MP_e.x_mp,MP_e.y_mp,self.bias_x,self.bias_y,self.Dh,
                                     self.Dh, self.dAx_grid_dt, self.dAy_grid_dt)
        #compute E-field in  boosted frame
        Ex_prime = dphi_prime_dx
        Ey_prime = dphi_prime_dy
        #compute E-field in  boosted frame
        Bx_prime = -1/(self.beta*c)*dAx_dt[0:MP_e.N_mp] - dAs_prime_dy
        By_prime = dAs_prime_dx + 1/(self.beta*c)*dAy_dt[0:MP_e.N_mp]
        Bz_prime = dAx_prime_dy - dAy_prime_dx

        #transform fields to lab frame
        Ex_sc_n = self.gamma*(Ex_prime + self.beta*c*By_prime)
        Ey_sc_n = self.gamma*(Ey_prime - self.beta*c*Bx_prime)

        Bx_sc_n = self.gamma*(Bx_prime - self.beta/c*Ey_prime)
        By_sc_n = self.gamma*(By_prime + self.beta/c*Ex_prime)
        Bz_sc_n = Bz_prime

        return Ex_sc_n, Ey_sc_n, Bx_sc_n, By_sc_n, Bz_sc_n

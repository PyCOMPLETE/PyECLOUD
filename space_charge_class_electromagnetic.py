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

    def __init__(self, chamb, Dh, gamma, Dt_sc=None, PyPICmode='FiniteDifferences_ShortleyWeller' , sparse_solver='scipy_slu',
                 f_telescope=None, target_grid=None, N_nodes_discard=None, N_min_Dh_main=None):

        super(space_charge_electromagnetic, self).__init__(chamb, Dh, Dt_sc, PyPICmode , sparse_solver,
                     f_telescope, target_grid, N_nodes_discard, N_min_Dh_main)

        self.flag_em_tracking = True

        self.state_Ax = self.PyPICobj.get_state_object()
        self.state_Ay = self.PyPICobj.get_state_object()

        self.state_Ax_old = None
        self.state_Ay_old = None

        # Ax and Ay at previous steps for computation of dAx_dz and dAy_dz
        self.Ax_old_grid = np.zeros((self.Nxg,self.Nyg))
        self.Ay_old_grid = np.zeros((self.Nxg,self.Nyg))
        self.dAx_grid_dt = np.zeros((self.Nxg,self.Nyg))
        self.dAy_grid_dt = np.zeros((self.Nxg,self.Nyg))

        self.gamma = gamma
        self.beta = np.sqrt(1-1/(gamma*gamma))

    def recompute_spchg_emfield(self, MP_e, flag_solve=True, flag_reset=True):
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

        #if not first passage compute derivatives
        if self.state_Ax_old != None and self.state_Ax_old != None:
            #compute time derivatives
            self.dAx_grid_dt = (self.state_Ax.phi -  self.state_Ax_old.phi)/self.Dt_sc
            self.dAy_grid_dt = (self.state_Ay.phi -  self.state_Ay_old.phi)/self.Dt_sc
        #if first passage set derivatives to zero
        else:
            self.dAx_grid_dt = np.zeros((self.Nxg,self.Nyg))
            self.dAy_grid_dt = np.zeros((self.Nxg,self.Nyg))

        self.state_Ax_old = self.state_Ax.get_state_object()
        self.state_Ay_old = self.state_Ay.get_state_object()

    def get_sc_em_field(self, MP_e):
        #compute un-primed potentials (with wrong sign)
        _, m_dAx_dy = self.state_Ax.gather(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp])
        m_dAy_dx, _ = self.state_Ay.gather(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp])
        m_dphi_dx, m_dphi_dy = self.PyPICobj.gather(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp])

        #fix signs and make primed
        dphi_prime_dx = -self.gamma*m_dphi_dx
        dphi_prime_dy = -self.gamma*m_dphi_dy
        dAx_prime_dy = -m_dAx_dy
        dAy_prime_dx = -m_dAy_dx
        #compute longitudinal magnetic potential
        dAs_prime_dx = -self.beta/c*dphi_prime_dx
        dAs_prime_dy = -self.beta/c*dphi_prime_dy
        # interpolate time derivatives to the particle positions
        dAx_dt, dAy_dt = iff.int_field(MP_e.x_mp,MP_e.y_mp,self.bias_x,self.bias_y,self.Dh,
                                     self.Dh, self.dAx_grid_dt, self.dAy_grid_dt)
        #compute E-field in  boosted frame
        Ex_prime = -dphi_prime_dx
        Ey_prime = -dphi_prime_dy
        #compute E-field in  boosted frame
        Bx_prime = dAs_prime_dy + 1/(self.beta*c)*dAy_dt[0:MP_e.N_mp]
        By_prime = -1/(self.beta*c)*dAx_dt[0:MP_e.N_mp] - dAs_prime_dx
        Bz_prime = dAy_prime_dx - dAx_prime_dy

        #transform fields to lab frame
        Ex_sc_n = self.gamma*(Ex_prime + self.beta*c*By_prime)
        Ey_sc_n = self.gamma*(Ey_prime - self.beta*c*Bx_prime)

        Bx_sc_n = self.gamma*(Bx_prime - self.beta/c*Ey_prime)
        By_sc_n = self.gamma*(By_prime + self.beta/c*Ex_prime)
        Bz_sc_n = Bz_prime

        #print(np.mean(np.abs(Ex_sc_n)/np.abs(MP_e.vy_mp[0:MP_e.N_mp]*Bz_sc_n-MP_e.vz_mp[0:MP_e.N_mp]*By_sc_n)))
        #print(np.mean(np.abs(Ex_sc_n)/(np.abs(c*Bz_sc_n-c*By_sc_n))))
        return Ex_sc_n, Ey_sc_n, Bx_sc_n, By_sc_n, Bz_sc_n

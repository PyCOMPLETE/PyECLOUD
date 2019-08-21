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
from scipy.constants import epsilon_0, mu_0
from scipy.constants import c  as c_light

import int_field_for as iff
import sys
from io import BytesIO as StringIO

import matplotlib.pyplot as plt
import matplotlib as mpl
import rhocompute as rhocom
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors
import int_field_for as iff

na = lambda x: np.array([x])

class space_charge_electromagnetic(space_charge, object):

    def __init__(self, chamb, Dh, gamma, Dt_sc=None, PyPICmode='FiniteDifferences_ShortleyWeller' , sparse_solver='scipy_slu',
                 f_telescope=None, target_grid=None, N_nodes_discard=None, N_min_Dh_main=None):

        super(space_charge_electromagnetic, self).__init__(chamb, Dh, Dt_sc, PyPICmode , sparse_solver,
                     f_telescope, target_grid, N_nodes_discard, N_min_Dh_main)

        self.flag_em_tracking = True
        # Initialize additional states for vector potential
        # The text trap makes sure that no output is shown in the terminal
        text_trap = StringIO()
        sys.stdout = text_trap
        self.state_Ax = self.PyPICobj.get_state_object()
        self.state_Ay = self.PyPICobj.get_state_object()
        sys.stdout = sys.__stdout__
        self.state_Ax_old = None
        self.state_Ay_old = None

        # Store the relativistic factors of the beam
        self.gamma = gamma
        self.beta = np.sqrt(1-1/(gamma*gamma))

    def recompute_spchg_emfield(self, MP_e, flag_solve=True, flag_reset=True):
        # Update the old states before scattering
        # The text trap makes sure that no output is shown in the terminal
        text_trap = StringIO()
        sys.stdout = text_trap
        self.state_Ax_old = self.state_Ax.get_state_object()
        self.state_Ay_old = self.state_Ay.get_state_object()
        sys.stdout = sys.__stdout__

        # Scatter the charge
        self.PyPICobj.scatter(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp], MP_e.nel_mp[0:MP_e.N_mp],
                charge=MP_e.charge, flag_add=not(flag_reset))

        # Scatter currents
        self.state_Ax.scatter(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp],
                epsilon_0 * mu_0 * MP_e.nel_mp[0:MP_e.N_mp] * MP_e.vx_mp[0:MP_e.N_mp],
                charge=MP_e.charge, flag_add=not(flag_reset))
        self.state_Ay.scatter(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp],
                epsilon_0 * mu_0 * MP_e.nel_mp[0:MP_e.N_mp] * MP_e.vy_mp[0:MP_e.N_mp],
                charge=MP_e.charge, flag_add=not(flag_reset))

        # Solve
        if flag_solve:
            self.PyPICobj.solve()
            self.PyPICobj.solve_states([self.state_Ax, self.state_Ay])

    def get_sc_em_field(self, MP_e):
        # Compute un-primed potentials (with wrong sign becase gather is meant to return E field..)
        _, m_dAx_dy = self.state_Ax.gather(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp])
        m_dAy_dx, _ = self.state_Ay.gather(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp])
        m_dphi_dx, m_dphi_dy = self.PyPICobj.gather(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp])
        # Fix signs
        dAx_dy = -m_dAx_dy
        dAy_dx = -m_dAy_dx

        # If not first passage compute time derivatives of Ax and Ay
        if self.state_Ax_old != None and self.state_Ax_old != None:
            dAx_dt = (self.state_Ax.gather_phi(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp]) -  self.state_Ax_old.gather_phi(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp]))/self.Dt_sc
            dAy_dt = (self.state_Ay.gather_phi(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp]) -  self.state_Ay_old.gather_phi(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp]))/self.Dt_sc
        # If first passage set derivatives to zero
        else:
            dAx_dt = np.zeros(MP_e.N_mp)
            dAy_dt = np.zeros(MP_e.N_mp)

        # Compute E-field in lab frame (E = -grad phi - dA/dt)
        Ex_sc_n = m_dphi_dx - dAx_dt
        Ey_sc_n = m_dphi_dy - dAy_dt

        # Compute B-field in lab frame (B = curl A)
        Bx_sc_n = 1/(self.beta*c_light)*dAy_dt
        By_sc_n = -1/(self.beta*c_light)*dAx_dt
        Bz_sc_n = dAy_dx - dAx_dy

        return Ex_sc_n, Ey_sc_n, Bx_sc_n, By_sc_n, Bz_sc_n

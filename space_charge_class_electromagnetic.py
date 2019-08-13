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

na = lambda x: np.array([x])

class space_charge_electromagnetic(space_charge, object):

    def __init__(self, chamb, Dh, gamma, Dt_sc=None, PyPICmode='FiniteDifferences_ShortleyWeller' , sparse_solver='scipy_slu',
                 f_telescope=None, target_grid=None, N_nodes_discard=None, N_min_Dh_main=None):

        super(space_charge_electromagnetic, self).__init__(chamb, Dh, Dt_sc, PyPICmode , sparse_solver,
                     f_telescope, target_grid, N_nodes_discard, N_min_Dh_main)

        self.flag_em_tracking = True
        text_trap = StringIO()
        sys.stdout = text_trap
        self.state_Ax = self.PyPICobj.get_state_object()
        self.state_Ay = self.PyPICobj.get_state_object()
        self.state_As = self.PyPICobj.get_state_object()
        sys.stdout = sys.__stdout__

        self.state_Ax_old = None
        self.state_Ay_old = None

        self.gamma = gamma
        self.beta = np.sqrt(1-1/(gamma*gamma))


        plt.figure(figsize=(14,7))

    def recompute_spchg_emfield(self, MP_e, flag_solve=True, flag_reset=True):
        #update the old states before scattering
        text_trap = StringIO()
        sys.stdout = text_trap
        self.state_Ax_old = self.state_Ax.get_state_object()
        self.state_Ay_old = self.state_Ay.get_state_object()
        sys.stdout = sys.__stdout__

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
            self.PyPICobj.solve_states([self.state_Ax, self.state_Ay, self.state_As])

    def get_sc_em_field(self, MP_e):
        #compute un-primed potentials (with wrong sign becase gather is meant to return E field..)
        _, m_dAx_dy = self.state_Ax.gather(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp])
        m_dAy_dx, _ = self.state_Ay.gather(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp])
        m_dphi_dx, m_dphi_dy = self.PyPICobj.gather(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp])
        #fix signs and make primed
        dphi_prime_dx = -self.gamma*m_dphi_dx
        dphi_prime_dy = -self.gamma*m_dphi_dy
        dAx_prime_dy = -m_dAx_dy
        dAy_prime_dx = -m_dAy_dx

        #compute longitudinal magnetic potential
        dAs_prime_dx = -self.beta/c_light*dphi_prime_dx
        dAs_prime_dy = -self.beta/c_light*dphi_prime_dy
        #compute E-field in  boosted frame
        Ex_prime = -dphi_prime_dx
        Ey_prime = -dphi_prime_dy

        #if not first passage compute derivatives of Ax and Ay
        if self.state_Ax_old != None and self.state_Ax_old != None:
            #compute time derivatives
            dAx_dt = (self.state_Ax.gather_phi(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp]) -  self.state_Ax_old.gather_phi(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp]))/self.Dt_sc
            dAy_dt = (self.state_Ay.gather_phi(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp]) -  self.state_Ay_old.gather_phi(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp]))/self.Dt_sc

        #if first passage set derivatives to zero
        else:
            dAx_dt = np.zeros(MP_e.N_mp)

            dAy_dt = np.zeros(MP_e.N_mp)

        dAx_prime_dt = dAx_dt
        dAy_prime_dt = dAy_dt
        #compute B-field in  boosted frame
        Bx_prime = dAs_prime_dy + 1/(self.gamma*self.beta*c_light)*dAy_prime_dt
        By_prime = -1/(self.gamma*self.beta*c_light)*dAx_prime_dt - dAs_prime_dx
        Bz_prime = dAy_prime_dx - dAx_prime_dy

        #transform fields to lab frame
        Ex_sc_n = self.gamma*(Ex_prime + self.beta*c_light*By_prime)
        Ey_sc_n = self.gamma*(Ey_prime - self.beta*c_light*Bx_prime)

        Bx_sc_n = self.gamma*(Bx_prime - self.beta/c_light*Ey_prime)
        By_sc_n = self.gamma*(By_prime + self.beta/c_light*Ex_prime)
        Bz_sc_n = Bz_prime

        ########################################################################
        dphi_dx = -m_dphi_dx
        dphi_dy = -m_dphi_dy
        self.compare_terms(dphi_dx, dphi_dy, dAx_dt, dAy_dt, Bx_sc_n, By_sc_n, self.beta, self.gamma, MP_e)
        ########################################################################

        return Ex_sc_n, Ey_sc_n, Bx_sc_n, By_sc_n, Bz_sc_n

    def compare_terms(self, dphi_dx_mp, dphi_dy_mp, dAx_dt_mp, dAy_dt_mp, Bx_sc_n_mp, By_sc_n_mp, beta, gamma, MP_e):
        xmin = np.min(self.xg)
        xmax = np.max(self.xg)
        ymin = np.min(self.yg)
        ymax = np.max(self.yg)
        dphi_dx = rhocom.compute_sc_rho(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp],dphi_dx_mp,self.bias_x,self.bias_y,self.Dh,self.Nxg,self.Nyg)
        dphi_dy = rhocom.compute_sc_rho(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp],dphi_dy_mp,self.bias_x,self.bias_y,self.Dh,self.Nxg,self.Nyg)
        dAx_dt = rhocom.compute_sc_rho(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp],dAx_dt_mp,self.bias_x,self.bias_y,self.Dh,self.Nxg,self.Nyg)
        dAy_dt = rhocom.compute_sc_rho(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp],dAy_dt_mp,self.bias_x,self.bias_y,self.Dh,self.Nxg,self.Nyg)
        Bx_sc_n = rhocom.compute_sc_rho(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp],Bx_sc_n_mp,self.bias_x,self.bias_y,self.Dh,self.Nxg,self.Nyg)
        By_sc_n = rhocom.compute_sc_rho(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp],By_sc_n_mp,self.bias_x,self.bias_y,self.Dh,self.Nxg,self.Nyg)
        #plt.plot(dphi_dx[inds],'b')
        fontsz = 15
        mpl.rc('xtick', labelsize=fontsz)
        mpl.rc('ytick', labelsize=fontsz)

        plt.subplot(2,3,3)
        plt.imshow(dphi_dx, cmap=None, norm=None, aspect='auto', interpolation=None,
                 alpha=None, origin='lower', extent=[xmin * 1e3, xmax * 1e3, ymin * 1e3, ymax * 1e3])
        plt.colorbar(format='%.0e')

        plt.title(r'$\frac{\partial \phi}{\partial x}$')
        plt.subplot(2,3,6)
        plt.imshow(dphi_dy, cmap=None, norm=None, aspect='auto', interpolation=None,
                 alpha=None, origin='lower', extent=[xmin * 1e3, xmax * 1e3, ymin * 1e3, ymax * 1e3])
        plt.title(r'$\frac{\partial \phi}{\partial y}$')
        plt.colorbar(format='%.0e')

        plt.subplot(2,3,1)
        ax = plt.gca()
        plt.imshow(dAx_dt, cmap=None, norm=None, aspect='auto', interpolation=None,
                alpha=None, origin='lower', extent=[xmin * 1e3, xmax * 1e3, ymin * 1e3, ymax * 1e3])
        plt.title(r'$\frac{\partial A_x}{\partial t}$')
        plt.colorbar(format='%.0e')

        plt.subplot(2,3,2)
        plt.imshow(-beta*c_light*By_sc_n, cmap=None, norm=None, aspect='auto', interpolation=None,
                 alpha=None, origin='lower', extent=[xmin * 1e3, xmax * 1e3, ymin * 1e3, ymax * 1e3])
        plt.title(r'$-\beta*c*B_y$')
        plt.colorbar(format='%.0e')

        plt.subplot(2,3,4)
        plt.imshow(dAy_dt, cmap=None, norm=None, aspect='auto', interpolation=None,
                 alpha=None, origin='lower', extent=[xmin * 1e3, xmax * 1e3, ymin * 1e3, ymax * 1e3])
        plt.title(r'$\frac{\partial A_y}{\partial t}$')
        plt.colorbar(format='%.0e')

        plt.subplot(2,3,5)
        plt.imshow(beta*c_light*Bx_sc_n, cmap=None, norm=None, aspect='auto', interpolation=None,
                 alpha=None, origin='lower', extent=[xmin * 1e3, xmax * 1e3, ymin * 1e3, ymax * 1e3])
        plt.title(r'$\beta*c*B_x$')
        plt.colorbar(format='%.0e')

        #plt.ylim(min(dAx_dt),max(dAx_dt))
        plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
        plt.draw()
        plt.pause(1e-5)
        plt.clf()
        #print('dAx_dt:%f'%(dAx_dt[10]))
        #print('-beta*c*By_sc_n:%f'%(-beta*c_light*By_sc_n[10]))
        #print('dphi_dx:%f'%(dphi_dx[10]))
        #print('err: %.16f'%(np.sqrt(np.mean(np.square(dAx_dt+beta*c_light*By_sc_n)))))

    '''
    def compare_terms(self, dphi_dx, dphi_dy, dAx_dt, dAy_dt, Bx_sc_n, By_sc_n, beta, gamma):
        indsx = sorted(range(len(dAx_dt)), key=lambda k: dAx_dt[k])
        indsy = sorted(range(len(dAy_dt)), key=lambda k: dAy_dt[k])
        indsphix = sorted(range(len(dphi_dx)), key=lambda k: dphi_dx[k])
        indsphiy = sorted(range(len(dphi_dy)), key=lambda k: dphi_dy[k])

        #plt.plot(dphi_dx[inds],'b')
        fontsz = 15
        mpl.rc('xtick', labelsize=fontsz)
        mpl.rc('ytick', labelsize=fontsz)
        plt.subplot(1,2,1)
        plt.plot(dAx_dt[indsx],'rx',label=r'$\frac{\partial A_x}{\partial t}$')
        plt.plot(-beta*c_light*By_sc_n[indsx],'b',label= r'$-\beta*c*B_y$')
        plt.plot(dAy_dt[indsy],'gx',label=r'$\frac{\partial A_y}{\partial t}$')
        plt.plot(beta*c_light*Bx_sc_n[indsy],'m',label= r'$\beta*c*B_x$')
        plt.xlabel('MP id', fontsize = fontsz)
        plt.ylabel(r'$[\frac{V}{m}]$', fontsize = fontsz)
        plt.legend(fontsize = fontsz)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), fontsize = fontsz)
        plt.subplot(1,2,2)
        plt.plot(dphi_dx[indsphix],'k',label= r'$\frac{\partial \phi}{\partial x}$')
        plt.plot(dphi_dy[indsphiy],'b',label= r'$\frac{\partial \phi}{\partial y}$')
        plt.xlabel('MP id', fontsize = fontsz)
        plt.ylabel(r'$[\frac{V}{m}]$', fontsize = fontsz)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), fontsize = fontsz)
        #plt.ylim(min(dAx_dt),max(dAx_dt))
        plt.legend(fontsize = fontsz)
        plt.draw()
        plt.pause(1e-5)
        plt.clf()
        #print('dAx_dt:%f'%(dAx_dt[10]))
        #print('-beta*c*By_sc_n:%f'%(-beta*c_light*By_sc_n[10]))
        #print('dphi_dx:%f'%(dphi_dx[10]))
        #print('err: %.16f'%(np.sqrt(np.mean(np.square(dAx_dt+beta*c_light*By_sc_n)))))
    '''

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
        self.ii = 0

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
        #fix signs
        dAx_dy = -m_dAx_dy
        dAy_dx = -m_dAy_dx

        #if not first passage compute derivatives of Ax and Ay
        if self.state_Ax_old != None and self.state_Ax_old != None:
            #compute time derivatives
            dAx_dt = (self.state_Ax.gather_phi(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp]) -  self.state_Ax_old.gather_phi(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp]))/self.Dt_sc
            dAy_dt = (self.state_Ay.gather_phi(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp]) -  self.state_Ay_old.gather_phi(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp]))/self.Dt_sc

        #if first passage set derivatives to zero
        else:
            dAx_dt = np.zeros(MP_e.N_mp)
            dAy_dt = np.zeros(MP_e.N_mp)

        #compute E-field in  boosted frame
        Ex_sc_n = m_dphi_dx - dAx_dt
        Ey_sc_n = m_dphi_dy - dAy_dt

        #compute B-field in  boosted frame
        Bx_sc_n = 1/(self.beta*c_light)*dAy_dt
        By_sc_n = -1/(self.beta*c_light)*dAx_dt
        Bz_sc_n = dAy_dx - dAx_dy

        if self.ii%20 == 0:
            #######################################################################
            dphi_dx = -m_dphi_dx
            dphi_dy = -m_dphi_dy
            self.compare_terms3(dphi_dx, dphi_dy, dAx_dt, dAy_dt, Bx_sc_n, By_sc_n, self.beta, self.gamma, MP_e)
            #######################################################################
        self.ii = self.ii+1
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
        fontsz = 20
        plt.set_cmap('bwr')

        plt.subplot(2,3,3)
        plt.imshow(dphi_dx.T, cmap=None, aspect='auto', interpolation=None,
                 alpha=None, origin='lower', extent=[xmin * 1e3, xmax * 1e3, ymin * 1e3, ymax * 1e3], norm=MidpointNormalize(midpoint=0))
        clb = plt.colorbar()
        clb.formatter.set_powerlimits((0, 0))
        clb.update_ticks()
        clb.set_label(r'$\left[\frac{V}{m}\right]$',fontsize = fontsz)
        plt.title(r'$\frac{\partial \phi}{\partial x}$',fontsize = fontsz)
        plt.subplot(2,3,6)
        plt.imshow(dphi_dy.T, cmap=None, norm=MidpointNormalize(midpoint=0), aspect='auto', interpolation=None,
                 alpha=None, origin='lower', extent=[xmin * 1e3, xmax * 1e3, ymin * 1e3, ymax * 1e3])
        plt.title(r'$\frac{\partial \phi}{\partial y}$',fontsize = fontsz)
        clb = plt.colorbar()
        clb.formatter.set_powerlimits((0, 0))
        clb.update_ticks()
        clb.set_label(r'$\left[\frac{V}{m}\right]$',fontsize = fontsz)
        plt.subplot(2,3,1)
        ax = plt.gca()
        plt.imshow(dAx_dt.T, cmap=None, aspect='auto', interpolation=None,
                alpha=None, origin='lower', extent=[xmin * 1e3, xmax * 1e3, ymin * 1e3, ymax * 1e3], norm=MidpointNormalize(midpoint=0))
        plt.title(r'$\frac{\partial A_x}{\partial t}$',fontsize = fontsz)
        clb = plt.colorbar()
        clb.formatter.set_powerlimits((0, 0))
        clb.update_ticks()
        clb.set_label(r'$\left[\frac{V}{m}\right]$',fontsize = fontsz)
        plt.subplot(2,3,2)
        plt.imshow(-beta*c_light*By_sc_n.T, cmap=None, norm=MidpointNormalize(midpoint=0), aspect='auto', interpolation=None,
                 alpha=None, origin='lower', extent=[xmin * 1e3, xmax * 1e3, ymin * 1e3, ymax * 1e3])
        plt.title(r'$-\beta*c*B_y$',fontsize = fontsz-5)
        clb = plt.colorbar()
        clb.formatter.set_powerlimits((0, 0))
        clb.update_ticks()
        clb.set_label(r'$\left[\frac{V}{m}\right]$',fontsize = fontsz)
        plt.subplot(2,3,4)
        plt.imshow(dAy_dt.T, cmap=None, norm=MidpointNormalize(midpoint=0), aspect='auto', interpolation=None,
                 alpha=None, origin='lower', extent=[xmin * 1e3, xmax * 1e3, ymin * 1e3, ymax * 1e3])
        plt.title(r'$\frac{\partial A_y}{\partial t}$',fontsize = fontsz)
        clb = plt.colorbar()
        clb.formatter.set_powerlimits((0, 0))
        clb.update_ticks()
        clb.set_label(r'$\left[\frac{V}{m}\right]$',fontsize = fontsz)
        plt.subplot(2,3,5)
        plt.imshow(beta*c_light*Bx_sc_n.T, cmap=None, norm=MidpointNormalize(midpoint=0), aspect='auto', interpolation=None,
                 alpha=None, origin='lower', extent=[xmin * 1e3, xmax * 1e3, ymin * 1e3, ymax * 1e3])
        plt.title(r'$\beta*c*B_x$',fontsize = fontsz-5)
        clb = plt.colorbar()
        clb.formatter.set_powerlimits((0, 0))
        clb.update_ticks()
        clb.set_label(r'$\left[\frac{V}{m}\right]$',fontsize = fontsz)
        #plt.ylim(min(dAx_dt),max(dAx_dt))
        plt.subplots_adjust(bottom=0.1, right=0.93, left = 0.05, top=0.9, hspace=0.5, wspace=0.5)
        #plt.draw()
        #plt.pause(1e-5)
        filename = str('comparison/%05d' % (self.ii)) + '.png'
        plt.savefig(filename, dpi=150)
        plt.clf()
        #print('dAx_dt:%f'%(dAx_dt[10]))
        #print('-beta*c*By_sc_n:%f'%(-beta*c_light*By_sc_n[10]))
        #print('dphi_dx:%f'%(dphi_dx[10]))
        #print('err: %.16f'%(np.sqrt(np.mean(np.square(dAx_dt+beta*c_light*By_sc_n)))))


    def compare_terms2(self, dphi_dx_mp, dphi_dy_mp, dAx_dt_mp, dAy_dt_mp, Bx_sc_n_mp, By_sc_n_mp, beta, gamma, MP_e):
        grad_phi_norm = np.sqrt(np.square(dphi_dx_mp)+np.square(dphi_dy_mp))
        dA_dt_norm = np.sqrt(np.square(dAx_dt_mp)+np.square(dAy_dt_mp))
        B_term_norm = beta*c_light*np.sqrt(np.square(Bx_sc_n_mp)+np.square(By_sc_n_mp))
        #plt.plot(dphi_dx[inds],'b')
        fontsz = 15
        mpl.rc('xtick', labelsize=fontsz)
        mpl.rc('ytick', labelsize=fontsz)

        fontsz = 20
        plt.set_cmap('bwr')

        win = 50

        plt.subplot(1,3,1)
        hist, bins = np.histogram(dA_dt_norm, bins = 1000)
        hist = np.true_divide(hist,sum(hist))
        hist = np.convolve(hist, np.ones((win,))/win, mode='same')
        plt.plot(bins[:-1],hist)
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.xlabel(r'$||\frac{d\mathbf{A}}{dt}||$')
        plt.subplot(1,3,2)
        hist, bins = np.histogram(B_term_norm, bins = 1000)
        hist = np.true_divide(hist,sum(hist))
        hist = np.convolve(hist, np.ones((win,))/win, mode='same')
        plt.plot(bins[:-1],hist)
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.xlabel(r'$||\beta c \hat{\mathbf{i}}_s \times \mathbf{A}||$')
        plt.subplot(1,3,3)
        hist, bins = np.histogram(grad_phi_norm, bins = 1000)
        hist = np.true_divide(hist,sum(hist))
        hist = np.convolve(hist, np.ones((win,))/win, mode='same')
        plt.plot(bins[:-1],hist)
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.xlabel(r'$||\nabla \phi||$')
        #plt.ylim(min(dAx_dt),max(dAx_dt))
        plt.subplots_adjust(bottom=0.1, right=0.93, left = 0.05, top=0.9, hspace=0.5, wspace=0.5)
        #plt.draw()
        #plt.pause(1e-5)
        filename = str('comparison_hists/%05d' % (self.ii)) + '.png'
        plt.savefig(filename, dpi=150)
        plt.clf()
        #print('dAx_dt:%f'%(dAx_dt[10]))
        #print('-beta*c*By_sc_n:%f'%(-beta*c_light*By_sc_n[10]))
        #print('dphi_dx:%f'%(dphi_dx[10]))
        #print('err: %.16f'%(np.sqrt(np.mean(np.square(dAx_dt+beta*c_light*By_sc_n)))))

    def compare_terms3(self, dphi_dx_mp, dphi_dy_mp, dAx_dt_mp, dAy_dt_mp, Bx_sc_n_mp, By_sc_n_mp, beta, gamma, MP_e):
        grad_phi_norm_mp = np.sqrt(np.square(dphi_dx_mp)+np.square(dphi_dy_mp))
        dA_dt_norm_mp = np.sqrt(np.square(dAx_dt_mp)+np.square(dAy_dt_mp))
        B_term_norm_mp = beta*c_light*np.sqrt(np.square(Bx_sc_n_mp)+np.square(By_sc_n_mp))
        scale=10000000
        grad_phi_norm_grid = rhocom.compute_sc_rho(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp],grad_phi_norm_mp,self.bias_x,self.bias_y,self.Dh,self.Nxg,self.Nyg)
        dA_dt_norm_grid = rhocom.compute_sc_rho(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp],dA_dt_norm_mp,self.bias_x,self.bias_y,self.Dh,self.Nxg,self.Nyg)
        B_term_norm_grid = rhocom.compute_sc_rho(MP_e.x_mp[0:MP_e.N_mp],MP_e.y_mp[0:MP_e.N_mp],B_term_norm_mp,self.bias_x,self.bias_y,self.Dh,self.Nxg,self.Nyg)
        sxx=0.0006656379575134914
        syy=0.0006849358075842791
        sxy=0
        #S = np.matrix([[sxx,sxy],[sxy,syy]])
        #m = [0,0]
        #XY_mp = np.random.multivariate_normal(m,S,1000)
        #X_mp = XY_mp[:,0]
        #Y_mp = XY_mp[:,1]
        #rho = np.linspace(0,3*sxx,1000)
        #theta = np.linspace(0,2*np.pi,1000)
        #X_mp = np.multiply(rho,np.cos(theta))
        #Y_mp = np.multiply(rho,np.sin(theta))
        xx = np.linspace(-3*sxx,3*sxx,100)
        yy = np.linspace(-3*syy,3*syy,100)
        X_mp = []
        Y_mp = []
        for i in range(len(xx)):
            for j in range(len(xx)):
                if xx[i]*xx[i]/(9*sxx*sxx)+yy[j]*yy[j]/(9*syy*syy)<=1:
                    X_mp.append(xx[i])
                    Y_mp.append(yy[j])

        grad_phi_norm = iff.int_field(X_mp,Y_mp,self.bias_x,self.bias_y,self.Dh,self.Dh,grad_phi_norm_grid,grad_phi_norm_grid)
        dA_dt_norm = iff.int_field(X_mp,Y_mp,self.bias_x,self.bias_y,self.Dh,self.Dh,dA_dt_norm_grid,dA_dt_norm_grid)
        B_term_norm = iff.int_field(X_mp,Y_mp,self.bias_x,self.bias_y,self.Dh,self.Dh,B_term_norm_grid,B_term_norm_grid)

        #plt.plot(dphi_dx[inds],'b')
        fontsz = 15
        mpl.rc('xtick', labelsize=fontsz)
        mpl.rc('ytick', labelsize=fontsz)

        fontsz = 20
        plt.set_cmap('bwr')

        win = 50

        #plt.subplot(1,3,1)
        hist, bins = np.histogram(dA_dt_norm, bins = 1000)
        hist = np.true_divide(hist,sum(hist))
        hist = np.convolve(hist, np.ones((win,))/win, mode='same')
        plt.semilogx(bins[:-1],hist,'b',label=r'$||\frac{d\mathbf{A}}{dt}||$')
        #plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        #plt.xlabel(r'$||\frac{d\mathbf{A}}{dt}||$', fontsize=fontsz)
        #plt.subplot(1,3,2)
        hist, bins = np.histogram(B_term_norm, bins = 1000)
        hist = np.true_divide(hist,sum(hist))
        hist = np.convolve(hist, np.ones((win,))/win, mode='same')
        plt.semilogx(bins[:-1],hist,'r', linestyle='--', dashes=(5, 10), label=r'$||\beta c \hat{\mathbf{i}}_s \times \mathbf{A}||$')
        #plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        #plt.xlabel(r'$||\beta c \hat{\mathbf{i}}_s \times \mathbf{A}||$', fontsize=fontsz)
        #plt.subplot(1,3,3)
        hist, bins = np.histogram(grad_phi_norm, bins = 1000)
        hist = np.true_divide(hist,sum(hist))
        hist = np.convolve(hist, np.ones((win,))/win, mode='same')
        plt.semilogx(bins[:-1],hist,'g',label=r'$||\nabla \phi||$')
        #plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        #plt.xlabel(r'$||\nabla \phi||$', fontsize=fontsz)
        #plt.ylim(min(dAx_dt),max(dAx_dt))
        plt.legend(fontsize=fontsz)
        plt.xlabel(r'$[N]$', fontsize=fontsz)
        #plt.subplots_adjust(bottom=0.15, right=0.93, left = 0.05, top=0.95, hspace=0.5, wspace=0.5)
        #plt.draw()
        #plt.pause(1e-5)
        filename = str('comparison_hists/%05d' % (self.ii)) + '.png'
        plt.savefig(filename, dpi=150)
        plt.clf()

    '''
    def compare_terms(self, dphi_dx, dphi_dy, dAx_dt, dAy_dt, Bx_sc_n, By_sc_n, beta, gamma):
        indsx = sorted(range(leIf False, the result will contain the number of samples in each bin. If True, the result is the value of the probability density function at the bin, normalized such that the integral over the range is 1. Note that the sum of the histogram values will not be equal to 1 unless bins of unity width are chosen; it is not a probability mass function.

n(dAx_dt)), key=lambda k: dAx_dt[k])
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

# set the colormap and centre the colorbar
class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

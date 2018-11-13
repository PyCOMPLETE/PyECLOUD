#-Begin-preamble-------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 7.6.0
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
from numpy.random import rand
import hist_for as histf
from scipy.constants import e, m_e


class MP_positions:
    def __init__(self, x, y, z):
        self.x_mp = x.copy()
        self.y_mp = y.copy()
        self.z_mp = z.copy()

class MP_system:
    def __init__(self, N_mp_max, nel_mp_ref_0, fact_split, fact_clean,
                 N_mp_regen_low, N_mp_regen, N_mp_after_regen,
                 Dx_hist_reg, Nx_reg, Ny_reg, Nvx_reg, Nvy_reg, Nvz_reg, regen_hist_cut, chamb,
                 N_mp_soft_regen=None, N_mp_after_soft_regen=None, charge=-e, mass=m_e):

        N_mp_max = int(N_mp_max)
        self.x_mp = np.zeros(N_mp_max, float)
        self.y_mp = np.zeros(N_mp_max, float)
        self.z_mp = np.zeros(N_mp_max, float)
        self.vx_mp = np.zeros(N_mp_max, float)
        self.vy_mp = np.zeros(N_mp_max, float)
        self.vz_mp = np.zeros(N_mp_max, float)
        self.nel_mp = np.zeros(N_mp_max, float)
        self.N_mp = 0

        self.nel_mp_ref = nel_mp_ref_0
        self.nel_mp_split = fact_split * self.nel_mp_ref
        #can be optimized (only true second)
        self.nel_mp_cl_th = fact_clean * self.nel_mp_ref
        self.nel_mp_ref_0 = nel_mp_ref_0
        self.fact_split = fact_split
        self.fact_clean = fact_clean

        self.N_mp_regen_low = N_mp_regen_low
        self.N_mp_regen = N_mp_regen
        self.N_mp_after_regen = N_mp_after_regen

        self.Nx_reg = Nx_reg
        self.Ny_reg = Ny_reg
        self.Nvx_reg = Nvx_reg
        self.Nvy_reg = Nvy_reg
        self.Nvz_reg = Nvz_reg
        self.chamb = chamb

        self.charge = charge
        self.mass = mass

        xg_hist_reg = np.arange(0, chamb.x_aper + 2. * Dx_hist_reg, Dx_hist_reg, float)
        xgr_hist_reg = xg_hist_reg[1:]
        xgr_hist_reg = xgr_hist_reg[::-1]#reverse array
        self.xg_hist_reg = np.concatenate((-xgr_hist_reg, xg_hist_reg), 0)
        self.Nxg_hist_reg = len(self.xg_hist_reg)
        self.bias_x_hist_reg = min(self.xg_hist_reg)
        self.Dx_hist_reg = Dx_hist_reg

        self.regen_hist_cut = regen_hist_cut

        self.flag_soft_regen = False
        if (N_mp_soft_regen is not None) and (N_mp_after_soft_regen is not None):
            self.flag_soft_regen = True
            self.N_mp_soft_regen = N_mp_soft_regen
            self.N_mp_after_soft_regen = N_mp_after_soft_regen

    def clean_small_MPs(self):

        print "Start clean. N_mp=%d Nel=%e"%(self.N_mp, np.sum(self.nel_mp[0:self.N_mp]))

        flag_clean = (self.nel_mp < self.nel_mp_cl_th)
        flag_keep = ~(flag_clean)
        flag_keep[self.N_mp:] = False
        self.N_mp = np.sum(flag_keep)

        self.x_mp[0:self.N_mp] = self.x_mp[flag_keep].copy()
        self.y_mp[0:self.N_mp] = self.y_mp[flag_keep].copy()
        self.z_mp[0:self.N_mp] = self.z_mp[flag_keep].copy()
        self.vx_mp[0:self.N_mp] = self.vx_mp[flag_keep].copy()
        self.vy_mp[0:self.N_mp] = self.vy_mp[flag_keep].copy()
        self.vz_mp[0:self.N_mp] = self.vz_mp[flag_keep].copy()
        self.nel_mp[0:self.N_mp] = self.nel_mp[flag_keep].copy()

        self.nel_mp[self.N_mp:] = 0.0

        print "Done clean. N_mp=%d Nel=%e"%(self.N_mp, np.sum(self.nel_mp[0:self.N_mp]))

    def set_nel_mp_ref(self, val):
        self.nel_mp_ref = val
        self.nel_mp_split = self.fact_split * val
        self.nel_mp_cl_th = self.fact_clean * val

    def check_for_soft_regeneration(self):

        if self.flag_soft_regen:

            if self.N_mp > self.N_mp_soft_regen:
                chrg = np.sum(self.nel_mp)
                erg = np.sum(0.5 / np.abs(self.charge / self.mass) * self.nel_mp[0:self.N_mp] * (self.vx_mp[0:self.N_mp] * self.vx_mp[0:self.N_mp] + self.vy_mp[0:self.N_mp] * self.vy_mp[0:self.N_mp] + self.vz_mp[0:self.N_mp] * self.vz_mp[0:self.N_mp]))

                new_nel_mp_ref = chrg / self.N_mp_after_soft_regen
                if new_nel_mp_ref < self.nel_mp_ref_0:
                    new_nel_mp_ref = self.nel_mp_ref_0

                #if new_nel_mp_ref>self.nel_mp_ref_0:removed from version 3.16
                print 'Start SOFT regeneration. N_mp=%d Nel_tot=%1.2e En_tot=%1.2e'%(self.N_mp, chrg, erg)

                self.set_nel_mp_ref(new_nel_mp_ref)

                death_prob = float(self.N_mp - self.N_mp_after_soft_regen) / float(self.N_mp)

                flag_keep = np.array(len(self.x_mp) * [False])
                flag_keep[:self.N_mp] = (rand(self.N_mp) > death_prob)
                self.N_mp = np.sum(flag_keep)

                self.x_mp[0:self.N_mp] = np.array(self.x_mp[flag_keep].copy())
                self.y_mp[0:self.N_mp] = np.array(self.y_mp[flag_keep].copy())
                self.z_mp[0:self.N_mp] = np.array(self.z_mp[flag_keep].copy())
                self.vx_mp[0:self.N_mp] = np.array(self.vx_mp[flag_keep].copy())
                self.vy_mp[0:self.N_mp] = np.array(self.vy_mp[flag_keep].copy())
                self.vz_mp[0:self.N_mp] = np.array(self.vz_mp[flag_keep].copy())
                self.nel_mp[0:self.N_mp] = np.array(self.nel_mp[flag_keep].copy())

                self.nel_mp[self.N_mp:] = 0.0

                chrg_before = chrg
                chrg_after = np.sum(self.nel_mp)

                correct_fact = chrg_before / chrg_after

                print 'Applied correction factor = %e'%correct_fact

                self.nel_mp[0:self.N_mp] = self.nel_mp[0:self.N_mp] * correct_fact

                chrg = np.sum(self.nel_mp)
                erg = np.sum(0.5 / np.abs(self.charge / self.mass) * self.nel_mp[0:self.N_mp] * (self.vx_mp[0:self.N_mp] * self.vx_mp[0:self.N_mp] + self.vy_mp[0:self.N_mp] * self.vy_mp[0:self.N_mp] + self.vz_mp[0:self.N_mp] * self.vz_mp[0:self.N_mp]))
                print 'Done SOFT regeneration. N_mp=%d Nel_tot=%1.2e En_tot=%1.2e'%(self.N_mp, chrg, erg)

    def check_for_regeneration(self):

        if (self.N_mp > self.N_mp_regen or (self.N_mp < self.N_mp_regen_low and self.nel_mp_ref > self.nel_mp_ref_0)):
            chrg = np.sum(self.nel_mp)
            erg = np.sum(0.5 / np.abs(self.charge / self.mass) * self.nel_mp[0:self.N_mp] * (self.vx_mp[0:self.N_mp] * self.vx_mp[0:self.N_mp] + self.vy_mp[0:self.N_mp] * self.vy_mp[0:self.N_mp] + self.vz_mp[0:self.N_mp] * self.vz_mp[0:self.N_mp]))
            print 'Start regeneration. N_mp=%d Nel_tot=%1.2e En_tot=%1.2e'%(self.N_mp, chrg, erg)

            new_nel_mp_ref = chrg / self.N_mp_after_regen
            if new_nel_mp_ref < self.nel_mp_ref_0:
                new_nel_mp_ref = self.nel_mp_ref_0

            self.set_nel_mp_ref(new_nel_mp_ref)

            hist_vect = np.zeros(self.Nxg_hist_reg, float)
            histf.compute_hist(self.x_mp[0:self.N_mp], self.nel_mp[0:self.N_mp], self.bias_x_hist_reg, self.Dx_hist_reg, hist_vect)
            nel_tot = np.sum(hist_vect)

            #eliminate the negligible part of the x histogram
            i_cut = 1
            cut_away = np.sum(hist_vect[:i_cut]) + np.sum(hist_vect[-i_cut:])
            while(cut_away < self.regen_hist_cut * nel_tot):
                i_cut = i_cut + 1
                cut_away = np.sum(hist_vect[:i_cut]) + np.sum(hist_vect[-i_cut:])

            x_max = (len(hist_vect) - i_cut + 1) * self.Dx_hist_reg + self.bias_x_hist_reg

            print 'x_max = %e'%x_max

            flag_clean = (abs(self.x_mp) > x_max)
            flag_keep = ~(flag_clean)
            flag_keep[self.N_mp:] = False
            self.N_mp = np.sum(flag_keep)

            self.x_mp[0:self.N_mp] = np.array(self.x_mp[flag_keep].copy())
            self.y_mp[0:self.N_mp] = np.array(self.y_mp[flag_keep].copy())
            self.z_mp[0:self.N_mp] = np.array(self.z_mp[flag_keep].copy())
            self.vx_mp[0:self.N_mp] = np.array(self.vx_mp[flag_keep].copy())
            self.vy_mp[0:self.N_mp] = np.array(self.vy_mp[flag_keep].copy())
            self.vz_mp[0:self.N_mp] = np.array(self.vz_mp[flag_keep].copy())
            self.nel_mp[0:self.N_mp] = np.array(self.nel_mp[flag_keep].copy())

            self.nel_mp[self.N_mp:] = 0.0

            #Assign particle to grid
            #
            vx_max = max(abs(self.vx_mp))
            vy_max = max(abs(self.vy_mp))
            vz_max = max(abs(self.vz_mp))
            #
            #
            Dx_reg = 2 * x_max / (self.Nx_reg - 1)
            bias_x = np.ceil(float(self.Nx_reg) / 2.)
            #Attention when trnslating to python

            Dy_reg = 2 * self.chamb.y_aper / (self.Ny_reg - 1)
            bias_y = np.ceil(float(self.Ny_reg) / 2.)
            #Attention when trnslating to python

            Dvx_reg = 2 * vx_max / (self.Nvx_reg - 1)
            bias_vx = np.ceil(float(self.Nvx_reg) / 2.)
            #Attention when trnslating to python

            Dvy_reg = 2 * vy_max / (self.Nvy_reg - 1)
            bias_vy = np.ceil(float(self.Nvy_reg) / 2.)
            #Attention when trnslating to python

            Dvz_reg = 2 * vz_max / (self.Nvz_reg - 1)
            bias_vz = np.ceil(float(self.Nvz_reg) / 2)
            #Attention when trnslating to python

            print 'particles_assigned_to grid'

            ##
            #% MATLAB-like indices
            ix_mp = np.around(self.x_mp[0:self.N_mp] / Dx_reg) + bias_x
            iy_mp = np.around(self.y_mp[0:self.N_mp] / Dy_reg) + bias_y
            ivx_mp = np.around(self.vx_mp[0:self.N_mp] / Dvx_reg) + bias_vx
            ivy_mp = np.around(self.vy_mp[0:self.N_mp] / Dvy_reg) + bias_vy
            ivz_mp = np.around(self.vz_mp[0:self.N_mp] / Dvz_reg) + bias_vz
            #
            #
            #
            indexes = (ix_mp - 1) * self.Ny_reg * self.Nvx_reg * self.Nvy_reg * self.Nvz_reg\
                          + (iy_mp - 1) * self.Nvx_reg * self.Nvy_reg * self.Nvz_reg\
                          + (ivx_mp - 1) * self.Nvy_reg * self.Nvz_reg\
                          + (ivy_mp - 1) * self.Nvz_reg\
                          + ivz_mp - 1
            indexes = np.int_(indexes)
            indices_nonzero_cells = np.array(list(set(indexes)))
            indices_nonzero_cells = np.sort(indices_nonzero_cells)

            vect_dens = dict(zip(indices_nonzero_cells, np.zeros(len(indices_nonzero_cells))))
            #lil_matrix((Nx_reg*Ny_reg*Nvx_reg*Nvy_reg*Nvz_reg,1));#allocate a sparse matrix
            #

            for i_mp in range(0, self.N_mp):
                index_curr = indexes[i_mp]
                vect_dens[index_curr] = vect_dens[index_curr] + self.nel_mp[i_mp]

            nonzero_cells = np.array(map(vect_dens.get, indices_nonzero_cells))

            #%% retrieve indices of nonempty cells
            #% NB use C-like indices
            divider = np.double(self.Ny_reg) * self.Nvx_reg * self.Nvy_reg * self.Nvz_reg
            ix_nonzero = np.int_(np.floor(indices_nonzero_cells / divider))
            indices_nonzero_cells = indices_nonzero_cells - ix_nonzero * divider

            divider = np.double(self.Nvx_reg) * self.Nvy_reg * self.Nvz_reg
            iy_nonzero = np.int_(np.floor(indices_nonzero_cells / divider))
            indices_nonzero_cells = indices_nonzero_cells - iy_nonzero * divider
            #
            divider = np.double(self.Nvy_reg) * self.Nvz_reg
            ivx_nonzero = np.int_(np.floor(indices_nonzero_cells / divider))
            indices_nonzero_cells = indices_nonzero_cells - ivx_nonzero * divider
            #
            ivy_nonzero = np.int_(np.floor(indices_nonzero_cells / np.double(self.Nvz_reg)))
            indices_nonzero_cells = indices_nonzero_cells - ivy_nonzero * self.Nvz_reg
            #
            ivz_nonzero = indices_nonzero_cells

            #
            #%pass to MATLAB-like indices
            ix_nonzero = ix_nonzero + 1
            iy_nonzero = iy_nonzero + 1
            ivx_nonzero = ivx_nonzero + 1
            ivy_nonzero = ivy_nonzero + 1
            ivz_nonzero = ivz_nonzero + 1

            #
            #
            x_nonzero = (ix_nonzero - bias_x) * Dx_reg
            y_nonzero = (iy_nonzero - bias_y) * Dy_reg
            vx_nonzero = (ivx_nonzero - bias_vx) * Dvx_reg
            vy_nonzero = (ivy_nonzero - bias_vy) * Dvy_reg
            vz_nonzero = (ivz_nonzero - bias_vz) * Dvz_reg
            #

            #%%

            num_MP_in_cell = nonzero_cells / self.nel_mp_ref

            intnum_MP_in_cell = np.int_(np.floor(num_MP_in_cell))
            rest = num_MP_in_cell - intnum_MP_in_cell
            ngen = len(rest)
            flag_rest = (rand(ngen) < rest)
            #
            intnum_MP_in_cell = intnum_MP_in_cell + np.int_(flag_rest)
            #% intnum_MP_in_cell_chk=intnum_MP_in_cell;
            N_mp_expect = np.sum(intnum_MP_in_cell)
            #
            #
            #
            self.x_mp = 0 * self.x_mp; self.y_mp = 0 * self.y_mp; self.z_mp = 0 * self.z_mp
            self.vx_mp = 0 * self.vx_mp; self.vy_mp = 0 * self.vy_mp; self.vz_mp = 0 * self.vz_mp
            self.nel_mp = 0 * self.nel_mp
            self.nel_mp[0:N_mp_expect] = np.ones(N_mp_expect) * self.nel_mp_ref
            #
            flag_add = (intnum_MP_in_cell > 0)
            n_add_step = np.sum(flag_add)
            #
            self.N_mp = 0
            while n_add_step > 0:

                x_temp = x_nonzero[flag_add] + Dx_reg * (rand(n_add_step) - 0.5)
                y_temp = y_nonzero[flag_add] + Dy_reg * (rand(n_add_step) - 0.5)

                vx_temp = vx_nonzero[flag_add] + Dvx_reg * (rand(n_add_step) - 0.5)
                vy_temp = vy_nonzero[flag_add] + Dvy_reg * (rand(n_add_step) - 0.5)
                vz_temp = vz_nonzero[flag_add] + Dvz_reg * (rand(n_add_step) - 0.5)

                x_nonzero_temp = x_nonzero[flag_add]
                y_nonzero_temp = y_nonzero[flag_add]

                flag_np = self.chamb.is_outside(x_temp, y_temp)
                #(((x_temp/x_aper)**2 + (y_temp/y_aper)**2)>=1)
                Nout = np.sum(flag_np)
                while(Nout > 0):
            #        indices=find(flag_np);
                    x_temp[flag_np] = x_nonzero_temp[flag_np] + Dx_reg * (rand(Nout) - 0.5)
                    y_temp[flag_np] = y_nonzero_temp[flag_np] + Dy_reg * (rand(Nout) - 0.5)
                    flag_np = self.chamb.is_outside(x_temp, y_temp)
                    #(((x_temp/x_aper)**2 + (y_temp/y_aper)**2)>=1)
                    Nout = np.sum(flag_np)

                self.x_mp[self.N_mp:self.N_mp + n_add_step] = x_temp
                self.y_mp[self.N_mp:self.N_mp + n_add_step] = y_temp

                self.vx_mp[self.N_mp:self.N_mp + n_add_step] = vx_temp
                self.vy_mp[self.N_mp:self.N_mp + n_add_step] = vy_temp
                self.vz_mp[self.N_mp:self.N_mp + n_add_step] = vz_temp
                self.N_mp = self.N_mp + n_add_step

                intnum_MP_in_cell[flag_add] = intnum_MP_in_cell[flag_add] - 1
                flag_add = intnum_MP_in_cell > 0
                n_add_step = np.sum(flag_add)

            #end

            chrg = np.sum(self.nel_mp)
            erg = np.sum(0.5 / np.abs(self.charge / self.mass) * self.nel_mp[0:self.N_mp] * (self.vx_mp[0:self.N_mp] * self.vx_mp[0:self.N_mp] + self.vy_mp[0:self.N_mp] * self.vy_mp[0:self.N_mp] + self.vz_mp[0:self.N_mp] * self.vz_mp[0:self.N_mp]))
            print 'Done regeneration. N_mp=%d Nel_tot=%1.2e En_tot=%1.2e'%(self.N_mp, chrg, erg)

    def add_uniform_MP_distrib(self, DNel, E_init, x_max, x_min, y_max, y_min):

            if x_max == None: x_max = self.chamb.x_aper
            if x_min == None: x_min = -self.chamb.x_aper
            if y_max == None: y_max = self.chamb.y_aper
            if y_min == None: y_min = -self.chamb.y_aper

            v0 = -np.sqrt(2. * (E_init / 3.) * np.abs(self.charge) / self.mass)

            N_new_MP = DNel / self.nel_mp_ref
            Nint_new_MP = int(np.floor(N_new_MP))
            rest = N_new_MP - Nint_new_MP
            Nint_new_MP = Nint_new_MP + int(rand() < rest)

            if Nint_new_MP > 0:

                x_temp = (x_max - x_min) * rand(Nint_new_MP) + x_min
                y_temp = (y_max - y_min) * rand(Nint_new_MP) + y_min

                flag_np = self.chamb.is_outside(x_temp, y_temp)#(((x_temp/x_aper)**2 + (y_temp/y_aper)**2)>=1);
                Nout = np.sum(flag_np)
                while(Nout > 0):
                    x_temp[flag_np] = (x_max - x_min) * rand(Nout) + x_min
                    y_temp[flag_np] = (y_max - y_min) * rand(Nout) + y_min
                    flag_np = self.chamb.is_outside(x_temp, y_temp)#(((x_temp/x_aper)**2 + (y_temp/y_aper)**2)>=1);
                    Nout = np.sum(flag_np)

                self.x_mp[self.N_mp:self.N_mp + Nint_new_MP] = x_temp
                #Be careful to the indexing when translating to python
                self.y_mp[self.N_mp:self.N_mp + Nint_new_MP] = y_temp
                self.z_mp[self.N_mp:self.N_mp + Nint_new_MP] = 0.
                #randn(Nint_new_MP,1)
                self.vx_mp[self.N_mp:self.N_mp + Nint_new_MP] = v0 * (rand() - 0.5)
                #if you note a towards down polarization look here
                self.vy_mp[self.N_mp:self.N_mp + Nint_new_MP] = v0 * (rand() - 0.5)
                self.vz_mp[self.N_mp:self.N_mp + Nint_new_MP] = v0 * (rand() - 0.5)
                self.nel_mp[self.N_mp:self.N_mp + Nint_new_MP] = self.nel_mp_ref

                self.N_mp = int(self.N_mp + Nint_new_MP)

    def add_uniform_ele_density(self, n_ele, E_init, x_max, x_min, y_max, y_min):

        if x_max is None:
            x_max = self.chamb.x_aper

        if x_min is None:
            x_min = -self.chamb.x_aper

        if y_max is None:
            y_max = self.chamb.y_aper

        if y_min is None:
            y_min = -self.chamb.y_aper

        v0 = -np.sqrt(2. * (E_init / 3.) * np.abs(self.charge) / self.mass)

        N_new_MP = n_ele * (x_max - x_min) * (y_max - y_min) / self.nel_mp_ref
        Nint_new_MP = int(np.floor(N_new_MP))
        rest = N_new_MP - Nint_new_MP
        Nint_new_MP = Nint_new_MP + int(rand() < rest)

        if Nint_new_MP > 0:

            x_temp = (x_max - x_min) * rand(Nint_new_MP) + x_min
            y_temp = (y_max - y_min) * rand(Nint_new_MP) + y_min

            flag_keep = ~self.chamb.is_outside(x_temp, y_temp)#(((x_temp/x_aper)**2 + (y_temp/y_aper)**2)>=1);
            x_temp = x_temp[flag_keep]
            y_temp = y_temp[flag_keep]
            Nint_new_MP = len(x_temp)

            self.x_mp[self.N_mp:self.N_mp + Nint_new_MP] = x_temp
            #Be careful to the indexing when translating to python
            self.y_mp[self.N_mp:self.N_mp + Nint_new_MP] = y_temp
            self.z_mp[self.N_mp:self.N_mp + Nint_new_MP] = 0.
            #randn(Nint_new_MP,1)
            self.vx_mp[self.N_mp:self.N_mp + Nint_new_MP] = v0 * (rand() - 0.5)
            #if you note a towards down polarization look here
            self.vy_mp[self.N_mp:self.N_mp + Nint_new_MP] = v0 * (rand() - 0.5)
            self.vz_mp[self.N_mp:self.N_mp + Nint_new_MP] = v0 * (rand() - 0.5)
            self.nel_mp[self.N_mp:self.N_mp + Nint_new_MP] = self.nel_mp_ref

            self.N_mp = int(self.N_mp + Nint_new_MP)

    def get_positions(self):
            return MP_positions(self.x_mp[:self.N_mp], self.y_mp[:self.N_mp], self.z_mp[:self.N_mp])

    def add_new_MPs(self, N_new_MP, nel_new_mp, x, y, z, vx, vy, vz):
        N_mp_old = self.N_mp
        N_mp_new = self.N_mp + N_new_MP
        self.x_mp[N_mp_old:N_mp_new] = x
        self.y_mp[N_mp_old:N_mp_new] = y
        self.z_mp[N_mp_old:N_mp_new] = z
        self.vx_mp[N_mp_old:N_mp_new] = vx
        self.vy_mp[N_mp_old:N_mp_new] = vy
        self.vz_mp[N_mp_old:N_mp_new] = vz
        self.nel_mp[N_mp_old:N_mp_new] = nel_new_mp
        self.N_mp = N_mp_new

    def add_from_file(self, filename_MPs):

        if type(filename_MPs) is str:
            import scipy.io as sio
            dict_MP_init = sio.loadmat(filename_MPs)
        else:
            dict_MP_init = filename_MPs

        Nint_new_MP = int(dict_MP_init['N_mp'])

        self.x_mp[self.N_mp:self.N_mp + Nint_new_MP] = np.squeeze(dict_MP_init['x_mp'])
        self.y_mp[self.N_mp:self.N_mp + Nint_new_MP] = np.squeeze(dict_MP_init['y_mp'])
        self.z_mp[self.N_mp:self.N_mp + Nint_new_MP] = np.squeeze(dict_MP_init['z_mp'])
        self.vx_mp[self.N_mp:self.N_mp + Nint_new_MP] = np.squeeze(dict_MP_init['vx_mp'])
        self.vy_mp[self.N_mp:self.N_mp + Nint_new_MP] = np.squeeze(dict_MP_init['vy_mp'])
        self.vz_mp[self.N_mp:self.N_mp + Nint_new_MP] = np.squeeze(dict_MP_init['vy_mp'])
        self.nel_mp[self.N_mp:self.N_mp + Nint_new_MP] = np.squeeze(dict_MP_init['nel_mp'])

        self.N_mp = int(self.N_mp + Nint_new_MP)

    def extract_dict(self):
        dict_MP = {
            'x_mp': self.x_mp[:self.N_mp].copy(),
            'y_mp': self.y_mp[:self.N_mp].copy(),
            'z_mp': self.z_mp[:self.N_mp].copy(),
            'vx_mp': self.vx_mp[:self.N_mp].copy(),
            'vy_mp': self.vy_mp[:self.N_mp].copy(),
            'vy_mp': self.vz_mp[:self.N_mp].copy(),
            'nel_mp': self.nel_mp[:self.N_mp].copy(),
            'N_mp': self.N_mp,
            }
        return dict_MP

    def init_from_dict(self, dict_MP):
        self.N_mp = 0
        self.add_from_file(dict_MP)


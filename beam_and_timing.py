#-Begin-preamble-------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 8.2.0
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
import scipy.io as sio
import scipy.special as sspe
from numpy import array
import int_field_for as iff


def bunch_train4(t, b_spac, t_offs, ppb_vect, sigmaz_vect):

    c = 299792458.
    zz = c * t
    N_bucket = len(ppb_vect)
    val = 0. * t

    for ii in range(0, N_bucket):
        if np.mod(ii, N_bucket / 20) == 0:
            print ('Beam profile generation %.0f'%(float(ii) / float(N_bucket) * 100) + """%""")

        ppb = ppb_vect[ii]
        sigmaz = sigmaz_vect[ii]
        if sigmaz > 0:
            z0 = c * (t_offs + ii * b_spac)
            mask_to_be_updated = (np.abs(zz - z0) < (10. * sigmaz))

            val[mask_to_be_updated] = val[mask_to_be_updated] + ppb / (sigmaz * np.sqrt(2 * np.pi)) *\
                np.exp(-(zz[mask_to_be_updated] - z0) * (zz[mask_to_be_updated] - z0) / (2 * sigmaz * sigmaz))

    return val


class beam_and_timing:
    def __init__(self, flag_bunched_beam, fact_beam, coast_dens, beam_charge, beam_field_file, lam_th_beam_field,
                 b_spac=None, sigmaz=None, t_offs=None, filling_pattern_file=None, Dt=None, t_end=None,
                 beam_long_prof_file=None, Dh_beam_field=None, f_telescope_beam=None, target_grid_beam=None,
                 N_nodes_discard_beam=None, N_min_Dh_main_beam=None,
                 chamb=None, sigmax=None, sigmay=None,
                 x_beam_pos=0., y_beam_pos=0., save_beam_field_file_as=None,
                 flag_secodary_beam=False, t_primary_beam=None,
                 Nx=None, Ny=None, nimag=None,
                 progress_mapgen_file=None):

        if chamb.is_outside(np.array([x_beam_pos]), np.array([y_beam_pos])):
            raise ValueError('The beam is outside the chamber!')

        flag_unif_Dt = True

        if flag_bunched_beam:
            print 'Start beam profile generation.'

            if type(filling_pattern_file) is str:
                if '.mat' in filling_pattern_file:
                    dict_fillp = sio.loadmat(filling_pattern_file)
                    ppb_vect = np.squeeze(dict_fillp['ppb_vect'])
                    if 'sigmaz_vect' in dict_fillp.keys():
                        sigmaz_vect = np.squeeze(dict_fillp['sigmaz_vect'])
                else:
                    raise ValueError('Input of filling scheme via txt files is discontinued!')
                    # Incompatible with Python3 because of the exec                   
                    # f = open(filling_pattern_file)
                    # exec(f.read())
                    # f.close()
            else:
                ppb_vect = np.atleast_1d(np.float_(np.array(filling_pattern_file)))

            if flag_secodary_beam:
                t = t_primary_beam
            else:
                N_bunches = len(ppb_vect)
                t = np.arange(0., N_bunches * b_spac + t_end + 2 * Dt, Dt)

            t_inter = t[-1]

            try:
                sigmaz_vect
            except NameError:
                if sigmaz != -1:
                    sigmaz_vect = 0 * ppb_vect
                    sigmaz_vect[ppb_vect > 0] = sigmaz
                else:
                    raise ValueError('Bunch length is not defined!!!')

            lam_t_array = bunch_train4(t, b_spac, t_offs, ppb_vect, sigmaz_vect)
            print 'Done beam profile generation.'
        else:
            print 'Loading beam profile from file:'
            print beam_long_prof_file

            dict_lam = sio.loadmat(beam_long_prof_file)
            t = np.squeeze(dict_lam['t'].real)
            lam_t_array = np.squeeze(dict_lam['lam_t_array'].real)
            t_inter = t[-1]

            Dt_vect = np.diff(t)
            if (np.max(np.abs(Dt_vect)) - np.mean(np.abs(Dt_vect))) > 1e-4 * np.mean(np.abs(Dt_vect)):
                flag_unif_Dt = False
            else:
                print 'Beam profile loaded from file.'
                print 'Uniform time step detected.'
                print 'The time step Dt provided in simulation_parameters.input will be ignored.'
                Dt = t[1] - t[0]
                print 'Time step set to Dt = %.3e s.'%Dt

        if flag_secodary_beam:
            if len(t) != len(t_primary_beam):
                raise ValueError('Time axes provided for primary and secondary beams should be identical!')
            if np.max(np.abs(t - t_primary_beam) / Dt) > 1e-4:
                raise ValueError('Time axes provided for primary and secondary beams should be identical!')

        Nt = len(t)
        lam_t_array = fact_beam * lam_t_array
        lam_t_array = lam_t_array + coast_dens

        N_pass_tot = int(np.ceil(t_inter / b_spac))

        flag_PyPIC_state_mode = False
        if beam_field_file == -1 or beam_field_file == 'computeFD':
            print 'No beam field file provided -> Calculate field using Poisson solver'

            if Dh_beam_field is None:
                raise ValueError('Grid size Dh_beam_field MUST be provided for beam field computation!')

            import space_charge_class as scc
            from numpy import exp, pi
            scb = scc.space_charge(chamb, Dh_beam_field, Dt_sc=1.)

            print 'Computing beam charge density'
            #rho=1./(2.*pi*sigmax*sigmay)*exp(-(scb.xn-x_beam_pos)**2/(2.*sigmax**2)-(scb.yn-y_beam_pos)**2/(2.*sigmay**2))

            rho = 1. / (4. * Dh_beam_field**2) * (sspe.erf((scb.xn - x_beam_pos) / (np.sqrt(2) * sigmax) + Dh_beam_field / (2 * np.sqrt(2) * sigmax))\
                                                  - sspe.erf((scb.xn - x_beam_pos) / (np.sqrt(2) * sigmax) - Dh_beam_field / (2 * np.sqrt(2) * sigmax)))\
                * (sspe.erf((scb.yn - y_beam_pos) / (np.sqrt(2) * sigmay) + Dh_beam_field / (2 * np.sqrt(2) * sigmay))\
                   - sspe.erf((scb.yn - y_beam_pos) / (np.sqrt(2) * sigmay) - Dh_beam_field / (2 * np.sqrt(2) * sigmay)))

            scb.compute_spchg_efield_from_rho(rho, flag_verbose=True)

            Ex_beam = scb.efx
            Ey_beam = scb.efy
            xx_beam = scb.xg
            yy_beam = scb.yg
            #self.scb = scb #######DEBUG!!
            print 'Done beam field computation.'

        elif beam_field_file == 'computeBE':

            print 'No beam field file provided -> Calculate field using Bassetti Erskine formula'

            if chamb.chamb_type != 'ellip':
                raise ValueError('You can only use Bassetti Erskine formula with an elliptic chamber!')
            if Nx is None or Ny is None or nimag is None:
                raise ValueError('Nx, Ny and nimag MUST be provided for Bassetti Erskine formula!')
            if x_beam_pos != 0. or y_beam_pos != 0.:
                raise ValueError('x_beam_pos, y_beam_pos and MUST be 0 for Bassetti Erskine formula!')

            print "sigmax=%.3e, sigmay=%.3e, Nx=%d, Ny=%d, nimag=%d"%(sigmax, sigmay, Nx, Ny, nimag)

            if progress_mapgen_file is not None:
                fprog = open(progress_mapgen_file, 'w')
                fprog.write('Bassetti Erskine\n')
                fprog.write("sigmax=%.3e, sigmay=%.3e, Nx=%d, Ny=%d, nimag=%d\n"%(sigmax, sigmay, Nx, Ny, nimag))
                fprog.close()

            import BassErsk as BE
            a = chamb.x_aper
            b = chamb.y_aper

            xmax = a * 1.02
            ymax = b * 1.02
            xx = np.linspace(-xmax, xmax, Nx)
            yy = np.linspace(-ymax, ymax, Ny)

            Ex = np.zeros((len(xx), len(yy)), dtype=complex)
            Ey = np.zeros((len(xx), len(yy)), dtype=complex)
            print 'Start beam field map generation.'
            for ii in range(len(xx)):

                if np.mod(ii, Nx / 20) == 0:
                    print ('Beam field map generation %.0f'%(float(ii) / float(Nx) * 100) + """%""")

                    if progress_mapgen_file is not None:
                        fprog = open(progress_mapgen_file, 'a')
                        fprog.write(('Done %.0f'%(float(ii) / float(Nx) * 100) + """%""" + '\n'))
                        fprog.close()

                for jj in range(len(yy)):
                    x = xx[ii]
                    y = yy[jj]
                    Ex_imag, Ey_imag = BE.ImageTerms(x, y, a, b, 0, 0, nimag)
                    Ex_BE, Ey_BE = BE.BassErsk(x, y, sigmax, sigmay)
                    Ex[ii, jj] = Ex_BE + Ex_imag
                    Ey[ii, jj] = Ey_BE + Ey_imag
                    Ex_beam = Ex.real
                    Ey_beam = Ey.real
                    xx_beam = xx
                    yy_beam = yy
                #to do: 1 - log progress somehow 2 - save file and check 3 - check secondary beams 4 - check the beam is centered
            if progress_mapgen_file is not None:
                fprog = open(progress_mapgen_file, 'a')
                fprog.write('Done.\n')
                fprog.close()
            print 'Done beam field map generation.'
        elif beam_field_file == 'compute_FDSW_multigrid':

            if Dh_beam_field is None:
                raise ValueError('Grid size Dh_beam_field MUST be provided for beam field computation!')
            if f_telescope_beam is None:
                raise ValueError('Aspect ratio MUST be provided for multigrid beam field computation!')
            if target_grid_beam is None:
                raise ValueError('Target grid MUST be provided for multigrid beam field computation!')
            if N_nodes_discard_beam is None:
                raise ValueError(' N_nodes_discard_beam MUST be provided for multigrid beam field computation!')
            if N_min_Dh_main_beam is None:
                raise ValueError(' N_min_Dh_main_beam MUST be provided for multigrid beam field computation!')

            import PyPIC.FiniteDifferences_ShortleyWeller_SquareGrid as PIC_FDSW
            PyPICmain = PIC_FDSW.FiniteDifferences_ShortleyWeller_SquareGrid(chamb=chamb, Dh=Dh_beam_field, sparse_solver='PyKLU')
            import PyPIC.MultiGrid as PIC_MG
            PyPICobj = PIC_MG.AddTelescopicGrids(pic_main=PyPICmain, f_telescope=f_telescope_beam, target_grid=target_grid_beam,
                                                 N_nodes_discard=N_nodes_discard_beam, N_min_Dh_main=N_min_Dh_main_beam, sparse_solver='PyKLU')

            # set rho
            # PyPICmain.rho=1./(2.*np.pi*sigmax*sigmay)*np.exp(-(PyPICmain.xn-x_beam_pos)**2/(2.*sigmax**2)-(PyPICmain.yn-y_beam_pos)**2/(2.*sigmay**2))
            PyPICmain.rho = 1. / (4. * Dh_beam_field**2) * (sspe.erf((PyPICmain.xn - x_beam_pos) / (np.sqrt(2) * sigmax) + Dh_beam_field / (2 * np.sqrt(2) * sigmax))\
                                                            - sspe.erf((PyPICmain.xn - x_beam_pos) / (np.sqrt(2) * sigmax) - Dh_beam_field / (2 * np.sqrt(2) * sigmax)))\
                * (sspe.erf((PyPICmain.yn - y_beam_pos) / (np.sqrt(2) * sigmay) + Dh_beam_field / (2 * np.sqrt(2) * sigmay))\
                   - sspe.erf((PyPICmain.yn - y_beam_pos) / (np.sqrt(2) * sigmay) - Dh_beam_field / (2 * np.sqrt(2) * sigmay)))
            for pic_dual in PyPICobj.pic_list:
                pic = pic_dual.pic_internal
                dh = pic_dual.pic_internal.Dh
                #pic.rho=1./(2.*np.pi*sigmax*sigmay)*np.exp(-(pic.xn-x_beam_pos)**2/(2.*sigmax**2)-(pic.yn-y_beam_pos)**2/(2.*sigmay**2))
                pic.rho = 1. / (4. * dh**2) * (sspe.erf((pic.xn - x_beam_pos) / (np.sqrt(2) * sigmax) + dh / (2 * np.sqrt(2) * sigmax))\
                                               - sspe.erf((pic.xn - x_beam_pos) / (np.sqrt(2) * sigmax) - dh / (2 * np.sqrt(2) * sigmax)))\
                    * (sspe.erf((pic.yn - y_beam_pos) / (np.sqrt(2) * sigmay) + dh / (2 * np.sqrt(2) * sigmay))\
                       - sspe.erf((pic.yn - y_beam_pos) / (np.sqrt(2) * sigmay) - dh / (2 * np.sqrt(2) * sigmay)))

            PyPICobj.solve()

            self.PyPIC_state = PyPICobj.get_state_object()
            self.get_beam_eletric_field = self._get_beam_eletric_field_PyPICstate

            flag_PyPIC_state_mode = True

            del(PyPICobj)

        else:
            print 'Loading beam field map from file:'
            print beam_field_file
            dict_beam = sio.loadmat(beam_field_file)

            Ex_beam = np.squeeze(dict_beam['Ex'].real)
            Ey_beam = np.squeeze(dict_beam['Ey'].real)
            xx_beam = np.squeeze(dict_beam['xx'].T)
            yy_beam = np.squeeze(dict_beam['yy'].T)

        if save_beam_field_file_as is not None:
            if flag_PyPIC_state_mode:
                raise ValueError('You cannot save the field maps in multigrid mode! Sorry...')
            sio.savemat(save_beam_field_file_as, {'xx': xx_beam, 'yy': yy_beam, 'Ex': Ex_beam, 'Ey': Ey_beam,\
                                                  'sigmax': sigmax, 'sigmay': sigmay,\
                                                  'x_aper': chamb.x_aper, 'y_aper': chamb.y_aper}, oned_as='row')

        if not flag_PyPIC_state_mode:

            xmin_beam = np.min(xx_beam)
            ymin_beam = np.min(yy_beam)
            dx_beam = xx_beam[1] - xx_beam[0]
            dy_beam = yy_beam[1] - yy_beam[0]

            self.Ex_beam = Ex_beam
            self.Ey_beam = Ey_beam
            self.xmin_beam = xmin_beam
            self.ymin_beam = ymin_beam
            self.dx_beam = dx_beam
            self.dy_beam = dy_beam
            self.xx_beam = xx_beam
            self.yy_beam = yy_beam

        if not(flag_unif_Dt):
            print 'Cloud simulation in non-uniform Dt mode.'
            print 'Dt provided in input will be used only as reference for savings and substeps.'

        self.flag_unif_Dt = flag_unif_Dt
        self.Nt = Nt
        self.b_spac = b_spac
        self.lam_t_array = lam_t_array
        self.beam_charge = beam_charge
        self.t = t
        self.Dt = Dt
        self.N_pass_tot = N_pass_tot

        self.lam_th_beam_field = lam_th_beam_field

        self.ii_curr = -1
        self.tt_curr = None
        self.lam_t_curr = None
        self.pass_numb = None
        self._pass_numb_old = -1

        self.sigmax = sigmax
        self.sigmay = sigmay
        self.x_beam_pos = x_beam_pos
        self.y_beam_pos = y_beam_pos

        self.flag_secodary_beam = flag_secodary_beam

    def next_time_step(self):
        self.ii_curr += 1
        self.tt_curr = self.t[self.ii_curr]
        self.Dt_curr = self.t[self.ii_curr + 1] - self.t[self.ii_curr]
        self.lam_t_curr = self.lam_t_array[self.ii_curr]
        self.pass_numb = int(np.floor(self.tt_curr / self.b_spac))
        self.flag_new_bunch_pass = (self.pass_numb > self._pass_numb_old)

        if self.flag_new_bunch_pass:
            self._pass_numb_old = self.pass_numb

    def end_simulation(self):
        return ((self.ii_curr + 2) >= self.Nt) # I need the last point to compute Dt

    def get_beam_eletric_field(self, MP_e):

        if (self.lam_t_curr > self.lam_th_beam_field) and (MP_e.N_mp > 0):
            ## compute beam electric field
            Ex_n_beam, Ey_n_beam = iff.int_field(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp],
                                                 self.xmin_beam, self.ymin_beam, self.dx_beam, self.dy_beam, self.Ex_beam, self.Ey_beam)
            Ex_n_beam = self.beam_charge * self.lam_t_curr * Ex_n_beam
            Ey_n_beam = self.beam_charge * self.lam_t_curr * Ey_n_beam

        else:
            Ex_n_beam = 0.
            Ey_n_beam = 0.

        return Ex_n_beam, Ey_n_beam

    def _get_beam_eletric_field_PyPICstate(self, MP_e):

        if (self.lam_t_curr > self.lam_th_beam_field) and (MP_e.N_mp > 0):
            ## compute beam electric field
            Ex_n_beam, Ey_n_beam = self.PyPIC_state.gather(MP_e.x_mp[0:MP_e.N_mp], MP_e.y_mp[0:MP_e.N_mp])
            Ex_n_beam = self.beam_charge * self.lam_t_curr * Ex_n_beam
            Ey_n_beam = self.beam_charge * self.lam_t_curr * Ey_n_beam
        else:
            Ex_n_beam = 0.
            Ey_n_beam = 0.

        return Ex_n_beam, Ey_n_beam

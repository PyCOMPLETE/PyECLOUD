#----------------------------------------------------------------------
#----------------------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 6.5.1
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
from scipy.constants import c, e, m_e

import myloadmat_to_obj as mlm

from geom_impact_ellip import ellip_cham_geom_object
from sec_emission_model_ECLOUD import SEY_model_ECLOUD
from sec_emission_model_accurate_low_ene import SEY_model_acc_low_ene
from sec_emission_model_ECLOUD_nunif import SEY_model_ECLOUD_non_unif
from sec_emission_model_cos_low_ener import SEY_model_cos_le
from sec_emission_model_flat_low_ener import SEY_model_flat_le
from sec_emission_model_from_file import SEY_model_from_file
import dynamics_Boris_f2py as dynB

import MP_system as MPs
import space_charge_class as scc
import impact_management_class as imc

import init

class MP_light(object):
    pass

class Ecloud(object):
    def __init__(self, L_ecloud, slicer, Dt_ref, pyecl_input_folder='./', flag_clean_slices=False,
                 slice_by_slice_mode=False, space_charge_obj=None, MP_e_mass=m_e, MP_e_charge=-e, **kwargs):


        print 'PyECLOUD Version 6.5.1'
        print 'PyHEADTAIL module'
        print 'Initializing ecloud from folder: '+pyecl_input_folder
        self.slicer = slicer
        self.Dt_ref = Dt_ref
        self.L_ecloud = L_ecloud

        self.pyecl_input_folder = pyecl_input_folder
        self.kwargs = kwargs

        config_dict = init.read_parameter_files(pyecl_input_folder, skip_beam_files=True)

        # Override config values with kwargs
        for attr, value in kwargs.items():
            if attr in ('x_beam_offset', 'y_beam_offset'):
                continue
            print('Ecloud init. From kwargs: %s = %r' % (attr, value))
            if attr in config_dict:
                config_dict[attr] = value
            else:
                print('Warning! What exactly does %s do? It is not part of any config file.' % attr)
                exec('%s=value') % attr

        cc = mlm.obj_from_dict(config_dict)

        #pyeclsaver=pysav.pyecloud_saver(logfile_path)


        if cc.switch_model=='ECLOUD_nunif':
            flag_non_unif_sey = 1
        else:
            flag_non_unif_sey = 0

        if cc.chamb_type=='ellip':
            chamb=ellip_cham_geom_object(cc.x_aper, cc.y_aper, flag_verbose_file=cc.flag_verbose_file)
        elif cc.chamb_type in ('polyg', 'polyg_cython'):
                import geom_impact_poly_fast_impact as gipfi
                chamb=gipfi.polyg_cham_geom_object(cc.filename_chm, flag_non_unif_sey, flag_verbose_file=cc.flag_verbose_file,
                                                   flag_verbose_stdout=cc.flag_verbose_stdout, flag_assume_convex=cc.flag_assume_convex)
        elif cc.chamb_type=='polyg_numpy':
            raise ValueError("chamb_type='polyg_numpy' not supported anymore")
            #~ chamb=gip.polyg_cham_geom_object(filename_chm, flag_non_unif_sey,
            #~ flag_verbose_file=flag_verbose_file, flag_verbose_stdout=flag_verbose_stdout)
        elif cc.chamb_type=='rect':
            import geom_impact_rect_fast_impact as girfi
            chamb = girfi.rect_cham_geom_object(cc.x_aper, cc.y_aper, flag_verbose_file=cc.flag_verbose_file,
                                                flag_verbose_stdout=cc.flag_verbose_stdout)
        else:
            raise ValueError('Chamber type not recognized (choose: ellip/rect/polyg)')


        MP_e=MPs.MP_system(cc.N_mp_max, cc.nel_mp_ref_0, cc.fact_split, cc.fact_clean,
                           cc.N_mp_regen_low, cc.N_mp_regen, cc.N_mp_after_regen,
                           cc.Dx_hist, cc.Nx_regen, cc.Ny_regen, cc.Nvx_regen, cc.Nvy_regen, cc.Nvz_regen, cc.regen_hist_cut, chamb,
                           N_mp_soft_regen=cc.N_mp_soft_regen, N_mp_after_soft_regen=cc.N_mp_after_soft_regen, charge=MP_e_charge,
                           mass=MP_e_mass)



        if cc.sparse_solver=='klu':
            print '''sparse_solver: 'klu' no longer supported --> going to PyKLU'''
            cc.sparse_solver='PyKLU'

        if space_charge_obj is not None:
            spacech_ele = space_charge_obj
        else:
            spacech_ele = scc.space_charge(chamb, cc.Dh_sc, Dt_sc=cc.Dt_sc, sparse_solver=cc.sparse_solver, PyPICmode=cc.PyPICmode,
                                           f_telescope=cc.f_telescope, target_grid=cc.target_grid, N_nodes_discard=cc.N_nodes_discard,
                                           N_min_Dh_main=cc.N_min_Dh_main)


        if cc.switch_model in (0, 'ECLOUD'):
            kwargs['flag_costheta_delta_scale'] = cc.flag_costheta_delta_scale
            kwargs['flag_costheta_Emax_shift'] = cc.flag_costheta_Emax_shift
            sey_mod = SEY_model_ECLOUD(cc.Emax, cc.del_max, cc.R0)
        elif cc.switch_model in (1, 'ACC_LOW'):
            sey_mod = SEY_model_acc_low_ene(cc.Emax, cc.del_max, cc.R0)
        elif cc.switch_model == 'ECLOUD_nunif':
            sey_mod = SEY_model_ECLOUD_non_unif(cc.chamb, cc.Emax, cc.del_max, cc.R0)
        elif cc.switch_model == 'cos_low_ene':
            sey_mod = SEY_model_cos_le(cc.Emax, cc.del_max, cc.R0)
        elif cc.switch_model == 'flat_low_ene':
                sey_mod = SEY_model_flat_le(cc.Emax, cc.del_max, cc.R0)
        elif cc.switch_model == 'perfect_absorber':
            sey_mod = None
        elif cc.switch_model == 'from_file':
            sey_mod = SEY_model_from_file(cc.sey_file, cc.flag_factor_costheta)


        flag_seg = (cc.flag_hist_impact_seg==1)

        if cc.switch_model=='perfect_absorber':
            import perfect_absorber_class as pac
            impact_man = pac.impact_management_perfect_absorber(
                cc.switch_no_increase_energy, chamb, sey_mod, cc.E_th, cc.sigmafit,cc.mufit, cc.Dx_hist, cc.scrub_en_th,
                cc.Nbin_En_hist, cc.En_hist_max, thresh_low_energy=cc.thresh_low_energy, flag_seg=flag_seg,
                cos_angle_width=cc.cos_angle_width, secondary_angle_distribution=cc.secondary_angle_distribution
            )
        else:
            impact_man = imc.impact_management(
                cc.switch_no_increase_energy, chamb, sey_mod, cc.E_th, cc.sigmafit, cc.mufit, cc.Dx_hist, cc.scrub_en_th,
                cc.Nbin_En_hist, cc.En_hist_max, thresh_low_energy=cc.thresh_low_energy, flag_seg=flag_seg,
                cos_angle_width=cc.cos_angle_width, secondary_angle_distribution=cc.secondary_angle_distribution
            )


        if cc.track_method == 'Boris':
            dynamics=dynB.pusher_Boris(cc.Dt, cc.B0x, cc.B0y, cc.B0z, cc.B_map_file, cc.fact_Bmap, cc.Bz_map_file, N_sub_steps=cc.N_sub_steps)
        #~ elif track_method == 'StrongBdip':
            #~ dynamics=dyndip.pusher_dipole_magnet(Dt,B)
        #~ elif track_method == 'StrongBgen':
            #~ dynamics=dyngen.pusher_strong_B_generalized(Dt, B0x, B0y,  \
                     #~ B_map_file, fact_Bmap, B_zero_thrhld)
        elif cc.track_method == 'BorisMultipole':
            import dynamics_Boris_multipole as dynmul
            dynamics=dynmul.pusher_Boris_multipole(Dt=cc.Dt, N_sub_steps=cc.N_sub_steps, B_multip=cc.B_multip)
        else:
            raise ValueError("""track_method should be 'Boris' or 'BorisMultipole' - others are not implemented in the PyEC4PyHT module""")


        if cc.init_unif_flag == 1:
            print("Adding inital %.2e electrons to the initial distribution") % cc.Nel_init_unif
            MP_e.add_uniform_MP_distrib(cc.Nel_init_unif, cc.E_init_unif, cc.x_max_init_unif, cc.x_min_init_unif,
                                        cc.y_max_init_unif, cc.y_min_init_unif)


        if cc.init_unif_edens_flag == 1:
            print("Adding inital %.2e electrons/m^3 to the initial distribution") % cc.init_unif_edens
            MP_e.add_uniform_ele_density(n_ele=cc.init_unif_edens, E_init=cc.E_init_unif_edens, x_max=cc.x_max_init_unif_edens,
                                         x_min=cc.x_min_init_unif_edens, y_max=cc.y_max_init_unif_edens, y_min=cc.y_min_init_unif_edens)


        if cc.filename_init_MP_state not in (-1, None):
            print("Adding inital electrons from: %s") % cc.filename_init_MP_state
            MP_e.add_from_file(cc.filename_init_MP_state)


        self.x_beam_offset = 0.
        self.y_beam_offset = 0.
        if 'x_beam_offset' in kwargs:
            self.x_beam_offset = kwargs['x_beam_offset']
        if 'y_beam_offset' in kwargs:
            self.y_beam_offset = kwargs['y_beam_offset']

        # initialize proton density probes
        self.save_ele_field_probes = False
        self.x_probes = -1
        self.y_probes = -1
        self.Ex_ele_last_track_at_probes = -1
        self.Ey_ele_last_track_at_probes = -1
        if 'probes_position' in kwargs.keys():
            self.save_ele_field_probes = True
            self.probes_position = kwargs['probes_position']
            self.N_probes = len(self.probes_position)
            self.x_probes = []
            self.y_probes = []
            for ii_probe in xrange(self.N_probes):
                self.x_probes.append(cc.probes_position[ii_probe]['x'])
                self.y_probes.append(cc.probes_position[ii_probe]['y'])

            self.x_probes = np.array(self.x_probes)
            self.y_probes = np.array(self.y_probes)

        self.N_tracks = 0

        spacech_ele.flag_decimate = False

        self.MP_e = MP_e
        self.dynamics = dynamics
        self.impact_man = impact_man
        self.spacech_ele = spacech_ele

        self.save_ele_distributions_last_track = False
        self.save_ele_potential_and_field = False
        self.save_ele_potential = False
        self.save_ele_field = False
        self.save_ele_MP_position = False
        self.save_ele_MP_velocity = False
        self.save_ele_MP_size = False

        self.track_only_first_time = False

        self.init_x = self.MP_e.x_mp[:self.MP_e.N_mp].copy()
        self.init_y = self.MP_e.y_mp[:self.MP_e.N_mp].copy()
        self.init_z = self.MP_e.z_mp[:self.MP_e.N_mp].copy()
        self.init_vx = self.MP_e.vx_mp[:self.MP_e.N_mp].copy()
        self.init_vy = self.MP_e.vy_mp[:self.MP_e.N_mp].copy()
        self.init_vz = self.MP_e.vz_mp[:self.MP_e.N_mp].copy()
        self.init_nel = self.MP_e.nel_mp[:self.MP_e.N_mp].copy()
        self.init_N_mp = self.MP_e.N_mp

        self.flag_clean_slices = flag_clean_slices

        self.slice_by_slice_mode = slice_by_slice_mode
        if self.slice_by_slice_mode:
            self.track = self._track_in_single_slice_mode
            self.finalize_and_reinitialize = self._finalize_and_reinitialize

    #    @profile
    def track(self, beam):

        if self.track_only_first_time:
            if self.N_tracks>0:
                print 'Warning: Track skipped because track_only_first_time is True.'
                return

        self._reinitialize()

        if hasattr(beam.particlenumber_per_mp, '__iter__'):
            raise ValueError('ecloud module assumes same size for all beam MPs')

        if self.flag_clean_slices:
            beam.clean_slices()

        slices = beam.get_slices(self.slicer)

        for i in xrange(slices.n_slices-1, -1, -1):

            # select particles in the slice
            ix = slices.particle_indices_of_slice(i)

            # slice size and time step
            dz = (slices.z_bins[i + 1] - slices.z_bins[i])

            self._track_single_slice(beam, ix, dz)

        self._finalize()

        self.N_tracks+=1

    def replace_with_recorded_field_map(self, delete_ecloud_data=True):

        if self.track_only_first_time:
            print 'Warning: replace_with_recorded_field_map resets track_only_first_time = False'
            self.track_only_first_time=False

        if not hasattr(self, 'efieldmap'):
            from Transverse_Efield_map_for_frozen_cloud import Transverse_Efield_map
            self.efieldmap = Transverse_Efield_map(xg = self.spacech_ele.xg, yg = self.spacech_ele.yg,
                Ex=self.Ex_ele_last_track, Ey=self.Ey_ele_last_track, L_interaction=self.L_ecloud,
                slicer = self.slicer,
                flag_clean_slices = True,
                x_beam_offset = self.x_beam_offset, y_beam_offset = self.y_beam_offset,
                slice_by_slice_mode = self.slice_by_slice_mode)

            self._ecloud_track = self.track

            self.track = self.efieldmap.track
            self.finalize_and_reinitialize = self.efieldmap.finalize_and_reinitialize

            if delete_ecloud_data:
                self.spacech_ele=None
                self.Mp_e = None
                self.init_nel = None
                self.init_vx = None
                self.init_vy = None
                self.init_vz = None
                self.init_x = None
                self.init_y = None
                self.init_z = None


        else:
            print 'Warning: efieldmap already exists. I do nothing.'

    def track_once_and_replace_with_recorded_field_map(self, bunch, delete_ecloud_data=True):
        self.save_ele_field = True
        self.track_only_first_time = True
        if self.slice_by_slice_mode:
            if not hasattr(bunch, '__iter__'):
                raise ValueError('A list of slices should be provided!')
            self._reinitialize()
            for slc in bunch:
                self.track(slc)
            self._finalize()
        else:
            self.track(bunch)
        self.save_ele_field = False
        self.track_only_first_time = False
        self.replace_with_recorded_field_map(delete_ecloud_data=delete_ecloud_data)


    def _track_single_slice(self, beam, ix, dz):

        MP_e = self.MP_e
        dynamics = self.dynamics
        impact_man = self.impact_man
        spacech_ele = self.spacech_ele

        dt = dz / (beam.beta * c)

        # define substep
        if dt>self.Dt_ref:
            N_sub_steps = int(np.round(dt/self.Dt_ref))
        else:
            N_sub_steps=1

        Dt_substep = dt/N_sub_steps
        #print Dt_substep, N_sub_steps, dt


        # beam field
        MP_p = MP_light()
        MP_p.x_mp = beam.x[ix]+self.x_beam_offset
        MP_p.y_mp = beam.y[ix]+self.y_beam_offset
        MP_p.nel_mp = beam.x[ix]*0.+beam.particlenumber_per_mp/dz#they have to become cylinders
        MP_p.N_mp = len(beam.x[ix])
        MP_p.charge = beam.charge
        #compute beam field (it assumes electrons!)
        spacech_ele.recompute_spchg_efield(MP_p)
        #scatter to electrons
        Ex_n_beam, Ey_n_beam = spacech_ele.get_sc_eletric_field(MP_e)


        ## compute electron field map
        spacech_ele.recompute_spchg_efield(MP_e)

        ## compute electron field on electrons
        Ex_sc_n, Ey_sc_n = spacech_ele.get_sc_eletric_field(MP_e)

        ## compute electron field on beam particles
        Ex_sc_p, Ey_sc_p = spacech_ele.get_sc_eletric_field(MP_p)

        ## Total electric field on electrons
        Ex_n=Ex_sc_n+Ex_n_beam;
        Ey_n=Ey_sc_n+Ey_n_beam;

        ## save position before motion step
        old_pos=MP_e.get_positions()

        ## motion electrons
        MP_e = dynamics.stepcustomDt(MP_e, Ex_n,Ey_n, Dt_substep=Dt_substep, N_sub_steps=N_sub_steps)

        ## impacts: backtracking and secondary emission
        MP_e = impact_man.backtrack_and_second_emiss(old_pos, MP_e)

        ## kick beam particles
        fact_kick = beam.charge/(beam.mass*beam.beta*beam.beta*beam.gamma*c*c)*self.L_ecloud
        beam.xp[ix]+=fact_kick*Ex_sc_p
        beam.yp[ix]+=fact_kick*Ey_sc_p

        if self.save_ele_distributions_last_track:
            self.rho_ele_last_track.append(spacech_ele.rho.copy())
            #print 'Here'

        if self.save_ele_potential:
            self.phi_ele_last_track.append(spacech_ele.phi.copy())

        if self.save_ele_field:
            self.Ex_ele_last_track.append(spacech_ele.efx.copy())
            self.Ey_ele_last_track.append(spacech_ele.efy.copy())

        if self.save_ele_MP_position:
            self.x_MP_last_track.append(MP_e.x_mp.copy())
            self.y_MP_last_track.append(MP_e.y_mp.copy())

        if self.save_ele_MP_velocity:
            self.vx_MP_last_track.append(MP_e.vx_mp.copy())
            self.vy_MP_last_track.append(MP_e.vy_mp.copy())

        if self.save_ele_MP_size:
            self.nel_MP_last_track.append(MP_e.nel_mp.copy())

        if self.save_ele_MP_position or self.save_ele_MP_velocity or self.save_ele_MP_size:
            self.N_MP_last_track.append(MP_e.N_mp)

        if self.save_ele_field_probes:
            MP_probes = MP_light()
            MP_probes.x_mp = self.x_probes
            MP_probes.y_mp = self.y_probes
            MP_probes.nel_mp = self.x_probes*0.+1. #fictitious charge of 1 C
            MP_probes.N_mp = len(self.x_probes)
            Ex_sc_probe, Ey_sc_probe = spacech_ele.get_sc_eletric_field(MP_probes)

            self.Ex_ele_last_track_at_probes.append(Ex_sc_probe.copy())
            self.Ey_ele_last_track_at_probes.append(Ey_sc_probe.copy())

    def _reinitialize(self):

        self.MP_e.x_mp[:self.init_N_mp] = self.init_x #it is a mutation and not a binding (and we have tested it :-))
        self.MP_e.y_mp[:self.init_N_mp] = self.init_y
        self.MP_e.z_mp[:self.init_N_mp] = self.init_z
        self.MP_e.vx_mp[:self.init_N_mp] = self.init_vx
        self.MP_e.vy_mp[:self.init_N_mp] = self.init_vy
        self.MP_e.vz_mp[:self.init_N_mp] = self.init_vz
        self.MP_e.nel_mp[:self.init_N_mp] = self.init_nel
        self.MP_e.N_mp = self.init_N_mp

        if self.save_ele_distributions_last_track:
            self.rho_ele_last_track = []

        if self.save_ele_potential_and_field:
            self.save_ele_potential = True
            self.save_ele_field = True

        if self.save_ele_potential:
            self.phi_ele_last_track = []

        if self.save_ele_field:
            self.Ex_ele_last_track = []
            self.Ey_ele_last_track = []

        if self.save_ele_MP_position:
            self.x_MP_last_track = []
            self.y_MP_last_track = []

        if self.save_ele_MP_velocity:
            self.vx_MP_last_track = []
            self.vy_MP_last_track = []

        if self.save_ele_MP_size:
            self.nel_MP_last_track = []

        if self.save_ele_MP_position or self.save_ele_MP_velocity or self.save_ele_MP_size:
            self.N_MP_last_track = []

        if self.save_ele_field_probes:
            self.Ex_ele_last_track_at_probes = []
            self.Ey_ele_last_track_at_probes = []

    def _finalize(self):

        if self.save_ele_distributions_last_track:
            self.rho_ele_last_track = np.array(self.rho_ele_last_track[::-1])

        if self.save_ele_potential:
            self.phi_ele_last_track = np.array(self.phi_ele_last_track[::-1])

        if self.save_ele_field:
            self.Ex_ele_last_track = np.array(self.Ex_ele_last_track[::-1])
            self.Ey_ele_last_track = np.array(self.Ey_ele_last_track[::-1])

        if self.save_ele_MP_position:
            self.x_MP_last_track = np.array(self.x_MP_last_track[::-1])
            self.y_MP_last_track = np.array(self.y_MP_last_track[::-1])

        if self.save_ele_MP_velocity:
            self.vx_MP_last_track = np.array(self.vx_MP_last_track[::-1])
            self.vy_MP_last_track = np.array(self.vy_MP_last_track[::-1])

        if self.save_ele_MP_size:
            self.nel_MP_last_track = np.array(self.nel_MP_last_track[::-1])

        if self.save_ele_MP_position or self.save_ele_MP_velocity or self.save_ele_MP_size:
                self.N_MP_last_track = np.array(self.N_MP_last_track[::-1])

        if self.save_ele_field_probes:
            self.Ex_ele_last_track_at_probes = np.array(self.Ex_ele_last_track_at_probes[::-1])
            self.Ey_ele_last_track_at_probes = np.array(self.Ey_ele_last_track_at_probes[::-1])

    def _finalize_and_reinitialize(self):
        self._finalize()
        self._reinitialize()

    def _track_in_single_slice_mode(self, beam):

        if hasattr(beam.particlenumber_per_mp, '__iter__'):
            raise ValueError('ecloud module assumes same size for all beam MPs')

        if self.flag_clean_slices:
            raise ValueError(
                    'track cannot clean the slices in slice-by-slice mode! ')

        if beam.slice_info is not 'unsliced':
            dz = beam.slice_info['z_bin_right']-beam.slice_info['z_bin_left']
            self._track_single_slice(beam, ix=np.arange(beam.macroparticlenumber), dz=dz)

    def generate_twin_ecloud_with_shared_space_charge(self):
        if hasattr(self, 'efieldmap'):
            raise ValueError('Ecloud has been replaced with field map. I cannot generate a twin ecloud!')
        return Ecloud(self.L_ecloud, self.slicer, self.Dt_ref, self.pyecl_input_folder, self.flag_clean_slices,
                self.slice_by_slice_mode, space_charge_obj=self.spacech_ele, **self.kwargs)


'''
def read_parameter_files_pyhdtl(pyecl_input_folder):
    switch_model=0
    simulation_param_file=pyecl_input_folder+'/simulation_parameters.input'

    save_mp_state_time_file = -1

    stopfile = 'stop'

    dec_fact_out = 1

    init_unif_flag = 0
    Nel_init_unif = None
    E_init_unif = 0.
    x_max_init_unif = None
    x_min_init_unif = None
    y_max_init_unif = None
    y_min_init_unif = None

    chamb_type = 'ellip'
    filename_chm = None

    x_aper = None
    y_aper = None
    flag_detailed_MP_info=0
    flag_hist_impact_seg = 0

    track_method= 'StrongBdip'

    B = 0.   #Tesla (if B=-1 computed from energy and bending radius)
    bm_totlen= -1 #m


    B0x = 0.
    B0y = 0.
    B0z = 0.
    B_map_file = None
    Bz_map_file = None
    N_sub_steps = 1
    fact_Bmap = 1.
    B_zero_thrhld = None


    # photoemission parameters
    photoem_flag = 0
    inv_CDF_refl_photoem_file = -1
    k_pe_st = -1
    refl_frac = -1
    alimit= -1
    e_pe_sigma = -1
    e_pe_max = -1
    x0_refl = -1
    y0_refl = -1
    out_radius = -1

    # gas ionization parameters
    gas_ion_flag = 0
    P_nTorr=-1
    sigma_ion_MBarn=-1
    Temp_K=-1
    unif_frac=-1
    E_init_ion=-1

    N_mp_soft_regen = None
    N_mp_after_soft_regen = None
    Dx = 0.
    Dy = 0.
    betafx = None
    betafy = None


    flag_verbose_file=False
    flag_verbose_stdout=False

    secondary_beams_file_list = []

    phem_resc_fac = 0.9999

    dec_fac_secbeam_prof=1

    el_density_probes=[]

    save_simulation_state_time_file = -1

    # detailed histogram
    x_min_hist_det=None
    x_max_hist_det=None
    y_min_hist_det=None
    y_max_hist_det=None
    Dx_hist_det=None

    sparse_solver = 'scipy_slu'

    B_multip = []

    PyPICmode = 'FiniteDifferences_ShortleyWeller'

    filename_init_MP_state = None

    # uniform initial density
    init_unif_edens_flag = 0
    init_unif_edens = None
    E_init_unif_edens= 0.
    x_max_init_unif_edens = None
    x_min_init_unif_edens = None
    y_max_init_unif_edens = None
    y_min_init_unif_edens = None

    flag_assume_convex = True


    # multigrid parameters
    f_telescope = None
    target_grid = None
    N_nodes_discard = None
    N_min_Dh_main = None


    f=open(simulation_param_file)
    exec(f.read())
    f.close()


    f=open(pyecl_input_folder+'/'+machine_param_file)
    exec(f.read())
    f.close()

    f=open(pyecl_input_folder+'/'+secondary_emission_parameters_file)
    exec(f.read())
    f.close()

    b_par = None# = pbf.beam_descr_from_fil(beam_parameters_file, betafx, Dx, betafy, Dy)

    flag_presence_sec_beams = False
    #~ if len(secondary_beams_file_list)>0:
        #~ flag_presence_sec_beams = True

    sec_b_par_list=[]
    #~ if flag_presence_sec_beams:
        #~ for sec_b_file in secondary_beams_file_list:
            #~ sec_b_par_list.append(pbf.beam_descr_from_fil(sec_b_file, betafx, Dx, betafy, Dy))

    if B==-1:
        B   = 2*pi*b_par.beta_rel*b_par.energy_J/(c*qe*bm_totlen)


    return b_par, x_aper, y_aper, B,\
    gas_ion_flag, P_nTorr, sigma_ion_MBarn, Temp_K, unif_frac, E_init_ion,\
    Emax, del_max, R0, E_th, sigmafit, mufit,\
    Dt, t_end, lam_th, t_ion, N_mp_max,\
    N_mp_regen, N_mp_after_regen, fact_split, fact_clean, nel_mp_ref_0,\
    Nx_regen, Ny_regen, Nvx_regen, Nvy_regen, Nvz_regen,regen_hist_cut,\
    N_mp_regen_low,\
    Dt_sc, Dh_sc, t_sc_ON,Dx_hist,r_center, scrub_en_th,\
    progress_path,  logfile_path, flag_movie, flag_sc_movie,\
    Dt_En_hist, Nbin_En_hist,En_hist_max, \
    photoem_flag, inv_CDF_refl_photoem_file, k_pe_st, refl_frac, alimit, e_pe_sigma,\
    e_pe_max,x0_refl, y0_refl, out_radius, \
    switch_model, switch_no_increase_energy, thresh_low_energy, save_mp_state_time_file, \
    init_unif_flag, Nel_init_unif, E_init_unif, x_max_init_unif, x_min_init_unif, y_max_init_unif, y_min_init_unif,\
    chamb_type, filename_chm, flag_detailed_MP_info, flag_hist_impact_seg,\
    track_method, B0x, B0y, B0z, B_map_file,  Bz_map_file, N_sub_steps, fact_Bmap, B_zero_thrhld,\
    N_mp_soft_regen, N_mp_after_soft_regen,\
    flag_verbose_file, flag_verbose_stdout,\
    flag_presence_sec_beams, sec_b_par_list, phem_resc_fac, dec_fac_secbeam_prof, el_density_probes, save_simulation_state_time_file,\
    x_min_hist_det, x_max_hist_det, y_min_hist_det, y_max_hist_det, Dx_hist_det, dec_fact_out, stopfile, sparse_solver, B_multip,\
    PyPICmode, filename_init_MP_state,\
    init_unif_edens_flag, init_unif_edens, E_init_unif_edens,\
    x_max_init_unif_edens, x_min_init_unif_edens, y_max_init_unif_edens, y_min_init_unif_edens, flag_assume_convex,\
    f_telescope, target_grid, N_nodes_discard, N_min_Dh_main
'''




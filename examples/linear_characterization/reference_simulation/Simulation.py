import os
# import Simulation_parameters as pp

import PyPARIS.communication_helpers as ch
import numpy as np
import PyPARIS.share_segments as shs
import time
import pickle
import h5py

from PyHEADTAIL.particles.slicing import UniformBinSlicer
from .sim_config_manager import SimConfig

class Simulation(object):

    def __init__(self, param_file='./Simulation_parameters.py'):

        self.pp = SimConfig(param_file)

    def init_all(self, generate_parent_eclouds=True,
            install_clouds=True):

        pp = self.pp

        self.N_turns = self.pp.N_turns

        self.n_slices = pp.n_slices

        # Prepare the machine without e-clouds
        self._build_machine()
        self._install_aperture()
        self._install_damper()
        self._install_impedance()

        # Split the machine
        self._split_machine_among_cores()

        # Generate and install e-clouds
        if generate_parent_eclouds:
            self._generate_parent_eclouds()
        if install_clouds:
            assert(generate_parent_eclouds)
            self._install_eclouds_in_machine_part()

        # Switch to footprint mode if needed
        if pp.footprint_mode:
            self._switch_to_footprint_mode()

    def init_master(self, generate_bunch=True, prepare_monitors=True):

        pp = self.pp

        # Manage multi-job operation
        if pp.footprint_mode:
            if pp.N_turns != pp.N_turns_target:
                raise ValueError(
                    "In footprint mode you need to set N_turns_target=N_turns_per_run!")
        self._setup_multijob_mode()


        # Define slicer
        self.slicer = UniformBinSlicer(
            n_slices=pp.n_slices, z_cuts=(-pp.z_cut, pp.z_cut)
        )

        # Prepare monitors
        if prepare_monitors:
            self._prepare_monitors()

        # generate the bunch and slice for the first turn
        if generate_bunch:
            self._generate_bunch()
            slice_obj_list = self.bunch.extract_slices(self.slicer)
            pieces_to_be_treated = slice_obj_list
        else:
            pieces_to_be_treated = []

        print("N_turns", self.N_turns)

        if pp.footprint_mode:
            self.recorded_particles = ParticleTrajectories(
                pp.n_macroparticles_for_footprint_track, self.N_turns
            )

        return pieces_to_be_treated

    def init_worker(self):
        pass

    def treat_piece(self, piece):
        for ele in self.mypart:
            ele.track(piece)

    def finalize_turn_on_master(self, pieces_treated):

        pp = self.pp

        # re-merge bunch
        self.bunch = sum(pieces_treated)

        # finalize present turn (with non parallel part, e.g. synchrotron motion)
        for ele in self.non_parallel_part:
            ele.track(self.bunch)

        # save results
        # print '%s Turn %d'%(time.strftime("%d/%m/%Y %H:%M:%S", time.localtime()), i_turn)
        self.bunch_monitor.dump(self.bunch)
        self.slice_monitor.dump(self.bunch)

        # prepare next turn (re-slice)
        new_pieces_to_be_treated = self.bunch.extract_slices(self.slicer)

        # order reset of all clouds
        orders_to_pass = ["reset_clouds"]

        # Save particles in case of footprint
        if pp.footprint_mode:
            self.recorded_particles.dump(self.bunch)

        # Check stop condition
        if self._check_stop_conditions():
            orders_to_pass.append("stop")
            self.SimSt.check_for_resubmit = False

        return orders_to_pass, new_pieces_to_be_treated

    def execute_orders_from_master(self, orders_from_master):
        if "reset_clouds" in orders_from_master:
            for ec in self.my_list_eclouds:
                ec.finalize_and_reinitialize()

    def finalize_simulation(self):

        pp = self.pp

        if pp.footprint_mode:
            # Get tunes
            from . import frequency_analysis as fa
            fa.get_tunes(self.recorded_particles,
                    filename_output='footprint.h5')

        else:
            # Finalize multijob info
            self._finalize_multijob_mode()

    def piece_to_buffer(self, piece):
        buf = ch.beam_2_buffer(piece)
        return buf

    def buffer_to_piece(self, buf):
        piece = ch.buffer_2_beam(buf)
        return piece


    def _build_machine(self):

        pp = self.pp

        self.optics_from_pickle = False

        if hasattr(pp, 'machine_class'):
            if pp.machine_class == 'Synchrotron':
                mode = 'synchrotron'
            elif pp.machine_class == 'LHC_custom':
                mode = 'LHC_custom'
            else:
                mode = 'custom_machine_class'
                raise ValueError('Not yet implemented')
        else:
            mode = 'LHC_custom'
            # kept as default for backward compatibility


        if mode == 'LHC_custom':

            # read the optics if needed
            if pp.optics_pickle_file is not None:
                with open(pp.optics_pickle_file) as fid:
                    optics = pickle.load(fid)
                    self.n_kick_smooth = np.sum(
                        ["_kick_smooth_" in nn for nn in optics["name"]]
                    )
                self.optics_from_pickle = True
            else:
                optics = None
                self.n_kick_smooth = pp.n_segments

            # define the machine
            from .LHC_custom import LHC

            self.machine = LHC(
                n_segments=pp.n_segments,
                machine_configuration=pp.machine_configuration,
                beta_x=pp.beta_x,
                beta_y=pp.beta_y,
                accQ_x=pp.Q_x,
                accQ_y=pp.Q_y,
                Qp_x=pp.Qp_x,
                Qp_y=pp.Qp_y,
                octupole_knob=pp.octupole_knob,
                optics_dict=optics,
                V_RF=pp.V_RF,
            )
        elif mode == 'synchrotron':
            from PyHEADTAIL.machines.synchrotron import Synchrotron
            self.machine = Synchrotron(
                optics_mode=pp.optics_mode,
                charge=pp.charge,
                mass=pp.mass,
                p0=pp.p0,
                circumference=pp.circumference,
                n_segments=pp.n_segments,
                name=pp.name,
                s=pp.s,
                alpha_x=pp.alpha_x,
                beta_x=pp.beta_x,
                D_x=pp.D_x,
                alpha_y=pp.alpha_y,
                beta_y=pp.beta_y,
                D_y=pp.D_y,
                accQ_x=pp.accQ_x,
                accQ_y=pp.accQ_y,
                Qp_x=pp.Qp_x,
                Qp_y=pp.Qp_y,
                app_x=pp.app_x,
                app_y=pp.app_y,
                app_xy=pp.app_xy,
                longitudinal_mode=pp.longitudinal_mode,
                Q_s=pp.Q_s,
                alpha_mom_compaction=pp.alpha_mom_compaction,
                h_RF=pp.h_RF,
                V_RF=pp.V_RF,
                dphi_RF=pp.dphi_RF,
                p_increment=pp.p_increment,
                RF_at=pp.RF_at,
                wrap_z=pp.wrap_z,
                other_detuners=pp.other_detuners,
            )

            if pp.optics_mode != 'smooth':
                raise ValueError('For arbitrary synchrotron only optics_mode="smooth" is implemented')

            self.n_kick_smooth = pp.n_segments
        else:
            raise ValueError('What?!')

        self.n_segments = self.machine.transverse_map.n_segments

        # compute sigma
        inj_opt = self.machine.transverse_map.get_injection_optics()
        sigma_x_inj = np.sqrt(inj_opt["beta_x"] * pp.epsn_x / self.machine.betagamma)
        sigma_y_inj = np.sqrt(inj_opt["beta_y"] * pp.epsn_y / self.machine.betagamma)

        if not self.optics_from_pickle:
            sigma_x_smooth = sigma_x_inj
            sigma_y_smooth = sigma_y_inj
        else:
            beta_x_smooth = None
            beta_y_smooth = None
            for ele in self.machine.one_turn_map:
                if ele in self.machine.transverse_map:
                    if "_kick_smooth_" in ele.name1:
                        if beta_x_smooth is None:
                            beta_x_smooth = ele.beta_x1
                            beta_y_smooth = ele.beta_y1
                        else:
                            if (
                                beta_x_smooth != ele.beta_x1
                                or beta_y_smooth != ele.beta_y1
                            ):
                                raise ValueError(
                                    "Smooth kicks must have all the same beta"
                                )

            if beta_x_smooth is None:
                sigma_x_smooth = None
                sigma_y_smooth = None
            else:
                sigma_x_smooth = np.sqrt(
                    beta_x_smooth * pp.epsn_x / self.machine.betagamma
                )
                sigma_y_smooth = np.sqrt(
                    beta_y_smooth * pp.epsn_y / self.machine.betagamma
                )

        self.sigma_x_inj = sigma_x_inj
        self.sigma_y_inj = sigma_y_inj
        self.sigma_x_smooth = sigma_x_smooth
        self.sigma_y_smooth = sigma_y_smooth

        self.n_non_parallelizable = 1 # longitudinal map


    def _generate_parent_eclouds(self):

        pp = self.pp

        sigma_x_inj = self.sigma_x_inj
        sigma_y_inj = self.sigma_y_inj
        sigma_x_smooth = self.sigma_x_smooth
        sigma_y_smooth = self.sigma_y_smooth

        # prepare e-cloud
        import PyECLOUD.PyEC4PyHT as PyEC4PyHT

        if pp.custom_target_grid_arcs is not None:
            target_grid_arcs = pp.custom_target_grid_arcs
        else:
            target_grid_arcs = {
                "x_min_target": -pp.target_size_internal_grid_sigma * sigma_x_smooth,
                "x_max_target": pp.target_size_internal_grid_sigma * sigma_x_smooth,
                "y_min_target": -pp.target_size_internal_grid_sigma * sigma_y_smooth,
                "y_max_target": pp.target_size_internal_grid_sigma * sigma_y_smooth,
                "Dh_target": pp.target_Dh_internal_grid_sigma * sigma_x_smooth,
            }
        self.target_grid_arcs = target_grid_arcs

        self.parent_eclouds = []

        nel_mp_ref_0 = (
                pp.init_unif_edens_dip * 4 * pp.x_aper * pp.y_aper
                / pp.N_MP_ele_init_dip
            )
        if pp.enable_arc_dip:
            # define MP size
            ecloud_dip = PyEC4PyHT.Ecloud(
                slice_by_slice_mode=True,
                L_ecloud=self.machine.circumference
                / self.n_kick_smooth
                * pp.fraction_device_dip,
                slicer=None,
                Dt_ref=pp.Dt_ref,
                pyecl_input_folder=pp.pyecl_input_folder,
                chamb_type=pp.chamb_type,
                x_aper=pp.x_aper,
                y_aper=pp.y_aper,
                filename_chm=pp.filename_chm,
                PyPICmode=pp.PyPICmode,
                Dh_sc=pp.Dh_sc_ext,
                N_min_Dh_main=pp.N_min_Dh_main,
                f_telescope=pp.f_telescope,
                N_nodes_discard=pp.N_nodes_discard,
                target_grid=target_grid_arcs,
                init_unif_edens_flag=pp.init_unif_edens_flag_dip,
                init_unif_edens=pp.init_unif_edens_dip,
                N_mp_max=pp.N_mp_max_dip,
                nel_mp_ref_0=nel_mp_ref_0,
                B_multip=pp.B_multip_dip,
                enable_kick_x=pp.enable_kick_x,
                enable_kick_y=pp.enable_kick_y,
                force_interp_at_substeps_interacting_slices=pp.force_interp_at_substeps_interacting_slices,
            )
            self.parent_eclouds.append(ecloud_dip)

        if pp.enable_arc_quad:
            ecloud_quad = PyEC4PyHT.Ecloud(
                slice_by_slice_mode=True,
                L_ecloud=self.machine.circumference
                / self.n_kick_smooth
                * pp.fraction_device_quad,
                slicer=None,
                Dt_ref=pp.Dt_ref,
                pyecl_input_folder=pp.pyecl_input_folder,
                chamb_type=pp.chamb_type,
                x_aper=pp.x_aper,
                y_aper=pp.y_aper,
                filename_chm=pp.filename_chm,
                PyPICmode=pp.PyPICmode,
                Dh_sc=pp.Dh_sc_ext,
                N_min_Dh_main=pp.N_min_Dh_main,
                f_telescope=pp.f_telescope,
                N_nodes_discard=pp.N_nodes_discard,
                target_grid=target_grid_arcs,
                N_mp_max=pp.N_mp_max_quad,
                nel_mp_ref_0=nel_mp_ref_0,
                B_multip=pp.B_multip_quad,
                filename_init_MP_state=pp.filename_init_MP_state_quad,
                enable_kick_x=pp.enable_kick_x,
                enable_kick_y=pp.enable_kick_y,
                force_interp_at_substeps_interacting_slices=pp.force_interp_at_substeps_interacting_slices,
            )
            self.parent_eclouds.append(ecloud_quad)

        if self.ring_of_CPUs.I_am_the_master and pp.enable_arc_dip:
            with open("multigrid_config_dip.txt", "w") as fid:
                if hasattr(ecloud_dip.spacech_ele.PyPICobj, "grids"):
                    fid.write(repr(ecloud_dip.spacech_ele.PyPICobj.grids))
                else:
                    fid.write("Single grid.")

            with open("multigrid_config_dip.pkl", "wb") as fid:
                if hasattr(ecloud_dip.spacech_ele.PyPICobj, "grids"):
                    pickle.dump(ecloud_dip.spacech_ele.PyPICobj.grids, fid)
                else:
                    pickle.dump("Single grid.", fid)

        if self.ring_of_CPUs.I_am_the_master and pp.enable_arc_quad:
            with open("multigrid_config_quad.txt", "w") as fid:
                if hasattr(ecloud_quad.spacech_ele.PyPICobj, "grids"):
                    fid.write(repr(ecloud_quad.spacech_ele.PyPICobj.grids))
                else:
                    fid.write("Single grid.")

            with open("multigrid_config_quad.pkl", "wb") as fid:
                if hasattr(ecloud_quad.spacech_ele.PyPICobj, "grids"):
                    pickle.dump(ecloud_quad.spacech_ele.PyPICobj.grids, fid)
                else:
                    pickle.dump("Single grid.", fid)

    def _install_damper(self):

        pp = self.pp

        if pp.enable_transverse_damper:
            # setup transverse damper
            from PyHEADTAIL.feedback.transverse_damper import TransverseDamper

            damper = TransverseDamper(
                dampingrate_x=pp.dampingrate_x, dampingrate_y=pp.dampingrate_y
            )
            self.machine.one_turn_map.append(damper)
            self.n_non_parallelizable += 1
            self.dampers = [damper]
        else:
            self.dampers = []

    def _install_aperture(self):

        pp = self.pp

        sigma_x_inj = self.sigma_x_inj
        sigma_y_inj = self.sigma_y_inj

        # setup transverse losses (to "protect" the ecloud)
        import PyHEADTAIL.aperture.aperture as aperture

        apt_xy = aperture.EllipticalApertureXY(
            x_aper=pp.target_size_internal_grid_sigma * sigma_x_inj,
            y_aper=pp.target_size_internal_grid_sigma * sigma_y_inj,
        )
        self.machine.one_turn_map.append(apt_xy)
        self.n_non_parallelizable += 1
        
        self.apertures = [apt_xy]


    def _split_machine_among_cores(self):

        pp = self.pp

        # We suppose that all the object that cannot
        # be slice parallelized are at the end of the ring
        i_end_parallel = len(self.machine.one_turn_map) - self.n_non_parallelizable

        # split the machine
        sharing = shs.ShareSegments(i_end_parallel, self.ring_of_CPUs.N_nodes)
        myid = self.ring_of_CPUs.myid
        i_start_part, i_end_part = sharing.my_part(myid)
        self.mypart = self.machine.one_turn_map[i_start_part:i_end_part]
        self.i_start_part = i_start_part
        if self.ring_of_CPUs.I_am_a_worker:
            print(
                "I am id=%d/%d (worker) and my part is %d long"
                % (myid, self.ring_of_CPUs.N_nodes, len(self.mypart))
            )
        elif self.ring_of_CPUs.I_am_the_master:
            self.non_parallel_part = self.machine.one_turn_map[i_end_parallel:]
            print(
                "I am id=%d/%d (master) and my part is %d long"
                % (myid, self.ring_of_CPUs.N_nodes, len(self.mypart))
            )

    def _install_eclouds_in_machine_part(self):
        # install eclouds in my part
        my_new_part = []
        self.my_list_eclouds = []
        for ele in self.mypart:
            my_new_part.append(ele)
            if ele in self.machine.transverse_map:
                if not self.optics_from_pickle or "_kick_smooth_" in ele.name1:
                    for ee in self.parent_eclouds:
                        ecloud_new = (
                            ee.generate_twin_ecloud_with_shared_space_charge()
                        )
                        my_new_part.append(ecloud_new)
                        self.my_list_eclouds.append(ecloud_new)
                elif (
                    "_kick_element_" in ele.name1 and pp.enable_eclouds_at_kick_elements
                ):

                    i_in_optics = list(optics["name"]).index(ele.name1)
                    kick_name = optics["name"][i_in_optics]
                    element_name = kick_name.split("_kick_element_")[-1]
                    L_curr = optics["L_interaction"][i_in_optics]

                    buildup_folder = pp.path_buildup_simulations_kick_elements.replace(
                        "!!!NAME!!!", element_name
                    )
                    chamber_fname = "%s_chamber.mat" % (element_name)

                    B_multip_curr = [0.0, optics["gradB"][i_in_optics]]

                    x_beam_offset = optics["x"][i_in_optics] * pp.orbit_factor
                    y_beam_offset = optics["y"][i_in_optics] * pp.orbit_factor

                    sigma_x_local = np.sqrt(
                        optics["beta_x"][i_in_optics]
                        * pp.epsn_x
                        / self.machine.betagamma
                    )
                    sigma_y_local = np.sqrt(
                        optics["beta_y"][i_in_optics]
                        * pp.epsn_y
                        / self.machine.betagamma
                    )

                    ecloud_ele = PyEC4PyHT.Ecloud(
                        slice_by_slice_mode=True,
                        L_ecloud=L_curr,
                        slicer=None,
                        Dt_ref=pp.Dt_ref,
                        pyecl_input_folder=pp.pyecl_input_folder,
                        chamb_type="polyg",
                        x_aper=None,
                        y_aper=None,
                        filename_chm=buildup_folder + "/" + chamber_fname,
                        PyPICmode=pp.PyPICmode,
                        Dh_sc=pp.Dh_sc_ext,
                        N_min_Dh_main=pp.N_min_Dh_main,
                        f_telescope=pp.f_telescope,
                        N_nodes_discard=pp.N_nodes_discard,
                        target_grid={
                            "x_min_target": -pp.target_size_internal_grid_sigma
                            * sigma_x_local
                            + x_beam_offset,
                            "x_max_target": pp.target_size_internal_grid_sigma
                            * sigma_x_local
                            + x_beam_offset,
                            "y_min_target": -pp.target_size_internal_grid_sigma
                            * sigma_y_local
                            + y_beam_offset,
                            "y_max_target": pp.target_size_internal_grid_sigma
                            * sigma_y_local
                            + y_beam_offset,
                            "Dh_target": pp.target_Dh_internal_grid_sigma
                            * sigma_y_local,
                        },
                        N_mp_max=pp.N_mp_max_quad,
                        nel_mp_ref_0=nel_mp_ref_0,
                        B_multip=B_multip_curr,
                        filename_init_MP_state=buildup_folder
                        + "/"
                        + pp.name_MP_state_file_kick_elements,
                        x_beam_offset=x_beam_offset,
                        y_beam_offset=y_beam_offset,
                        enable_kick_x=pp.enable_kick_x,
                        enable_kick_y=pp.enable_kick_y,
                        force_interp_at_substeps_interacting_slices=pp.force_interp_at_substeps_interacting_slices,
                    )

                    my_new_part.append(ecloud_ele)
                    self.my_list_eclouds.append(ecloud_ele)

        self.mypart = my_new_part

    def _install_impedance(self):

        pp = self.pp
        if hasattr(pp, 'enable_impedance'):
            if pp.enable_impedance:

                slicer_for_wakefields = UniformBinSlicer(
                        pp.n_slices_wake, z_cuts=(-pp.z_cut, pp.z_cut))

                import PyHEADTAIL.impedances.wakes as wakes
                wake = wakes.CircularResonator(R_shunt=pp.resonator_R_shunt,
                        frequency=pp.resonator_frequency,
                        Q=pp.resonator_Q)
                wake_element = wakes.WakeField(slicer_for_wakefields, wake)
                self.machine.one_turn_map.append(wake_element)
                self.n_non_parallelizable += 1
                self.impedances = [wake_element]
            else:
                self.impedances = []

    def _switch_to_footprint_mode(self):

        pp = self.pp

        print("Proc. %d computing maps" % self.ring_of_CPUs.myid)
        # generate a bunch
        bunch_for_map = self.machine.generate_6D_Gaussian_bunch_matched(
            n_macroparticles=pp.n_macroparticles_for_footprint_map,
            intensity=pp.intensity,
            epsn_x=pp.epsn_x,
            epsn_y=pp.epsn_y,
            sigma_z=pp.sigma_z,
        )

        # Slice the bunch
        slicer_for_map = UniformBinSlicer(
            n_slices=pp.n_slices, z_cuts=(-pp.z_cut, pp.z_cut)
        )
        slices_list_for_map = bunch_for_map.extract_slices(slicer_for_map)

        # Track the previous part of the machine
        for ele in self.machine.one_turn_map[:self.i_start_part]:
            for ss in slices_list_for_map:
                ele.track(ss)

        # Measure optics, track and replace clouds with maps
        list_ele_type = []
        list_meas_beta_x = []
        list_meas_alpha_x = []
        list_meas_beta_y = []
        list_meas_alpha_y = []
        for ele in self.mypart:
            list_ele_type.append(str(type(ele)))
            # Measure optics
            bbb = sum(slices_list_for_map)
            list_meas_beta_x.append(bbb.beta_Twiss_x())
            list_meas_alpha_x.append(bbb.alpha_Twiss_x())
            list_meas_beta_y.append(bbb.beta_Twiss_y())
            list_meas_alpha_y.append(bbb.alpha_Twiss_y())

            if ele in self.my_list_eclouds:
                ele.track_once_and_replace_with_recorded_field_map(
                    slices_list_for_map
                )
            else:
                for ss in slices_list_for_map:
                    ele.track(ss)
        print("Proc. %d done with maps" % self.ring_of_CPUs.myid)

        with open("measured_optics_%d.pkl" % self.ring_of_CPUs.myid, "wb") as fid:
            pickle.dump(
                {
                    "ele_type": list_ele_type,
                    "beta_x": list_meas_beta_x,
                    "alpha_x": list_meas_alpha_x,
                    "beta_y": list_meas_beta_y,
                    "alpha_y": list_meas_alpha_y,
                },
                fid,
            )

        # remove RF
        if self.ring_of_CPUs.I_am_the_master:
            self.non_parallel_part.remove(self.machine.longitudinal_map)

    def _generate_bunch(self):

        pp = self.pp

        # generate a bunch
        if pp.footprint_mode:
            self.bunch = self.machine.generate_6D_Gaussian_bunch_matched(
                n_macroparticles=pp.n_macroparticles_for_footprint_track,
                intensity=pp.intensity,
                epsn_x=pp.epsn_x,
                epsn_y=pp.epsn_y,
                sigma_z=pp.sigma_z,
            )
        elif self.SimSt.first_run:

            if pp.bunch_from_file is not None:
                print("Loading bunch from file %s ..." % pp.bunch_from_file)
                with h5py.File(pp.bunch_from_file, "r") as fid:
                    self.bunch = self.buffer_to_piece(np.array(fid["bunch"]).copy())
                print("Bunch loaded from file.\n")

            else:
                self.bunch = self.machine.generate_6D_Gaussian_bunch_matched(
                    n_macroparticles=pp.n_macroparticles,
                    intensity=pp.intensity,
                    epsn_x=pp.epsn_x,
                    epsn_y=pp.epsn_y,
                    sigma_z=pp.sigma_z,
                )

                # Recenter all slices
                if hasattr(pp, 'recenter_all_slices'):
                    if pp.recenter_all_slices:
                        print('Recentering all slices')
                        temp_slices = self.bunch.get_slices(self.slicer)
                        for ii in range(temp_slices.n_slices):
                            ix = temp_slices.particle_indices_of_slice(ii)
                            if len(ix) > 0:
                                self.bunch.x[ix] -= np.mean(self.bunch.x[ix])
                                self.bunch.xp[ix] -= np.mean(self.bunch.xp[ix])
                                self.bunch.y[ix] -= np.mean(self.bunch.y[ix])
                                self.bunch.yp[ix] -= np.mean(self.bunch.yp[ix])

                # compute initial displacements
                inj_opt = self.machine.transverse_map.get_injection_optics()
                sigma_x = np.sqrt(
                    inj_opt["beta_x"] * pp.epsn_x / self.machine.betagamma
                )
                sigma_y = np.sqrt(
                    inj_opt["beta_y"] * pp.epsn_y / self.machine.betagamma
                )
                x_kick = pp.x_kick_in_sigmas * sigma_x
                y_kick = pp.y_kick_in_sigmas * sigma_y

                # apply initial displacement
                if not pp.footprint_mode:
                    self.bunch.x += x_kick
                    self.bunch.y += y_kick

                print("Bunch initialized.")
        else:
            print("Loading bunch from file...")
            with h5py.File(
                "bunch_status_part%02d.h5" % (self.SimSt.present_simulation_part - 1), "r"
            ) as fid:
                self.bunch = self.buffer_to_piece(np.array(fid["bunch"]).copy())
            print("Bunch loaded from file.")

    def _prepare_monitors(self):

        pp = self.pp

        if hasattr(pp, 'write_buffer_every'):
            write_buffer_every = pp.write_buffer_every
        else:
            write_buffer_every = 3

        # define a bunch monitor
        from PyHEADTAIL.monitors.monitors import BunchMonitor

        self.bunch_monitor = BunchMonitor(
            "bunch_evolution_%02d" % self.SimSt.present_simulation_part,
            pp.N_turns,
            {"Comment": "PyHDTL simulation"},
            write_buffer_every=write_buffer_every,
        )

        # define a slice monitor
        from PyHEADTAIL.monitors.monitors import SliceMonitor

        self.slice_monitor = SliceMonitor(
            "slice_evolution_%02d" % self.SimSt.present_simulation_part,
            pp.N_turns,
            self.slicer,
            {"Comment": "PyHDTL simulation"},
            write_buffer_every=write_buffer_every,
        )

    def _setup_multijob_mode(self):

        pp = self.pp

        check_for_resubmit = True
        if hasattr(pp, "check_for_resubmit"):
            check_for_resubmit = pp.check_for_resubmit
        import PyPARIS_sim_class.Save_Load_Status as SLS

        SimSt = SLS.SimulationStatus(
            N_turns_per_run=pp.N_turns,
            check_for_resubmit=check_for_resubmit,
            N_turns_target=pp.N_turns_target,
        )
        SimSt.before_simulation()
        self.SimSt = SimSt

    def _check_stop_conditions(self):

        pp = self.pp

        stop = False
        # check if simulation has to be stopped
        # 1. for beam losses
        if (
            not pp.footprint_mode
            and self.bunch.macroparticlenumber < pp.sim_stop_frac * pp.n_macroparticles
        ):
            stop = True
            print("Stop simulation due to beam losses.")

        # 2. for the emittance growth
        if pp.flag_check_emittance_growth:
            epsn_x_max = (pp.epsn_x) * (1 + pp.epsn_x_max_growth_fraction)
            epsn_y_max = (pp.epsn_y) * (1 + pp.epsn_y_max_growth_fraction)
            if not pp.footprint_mode and (
                self.bunch.epsn_x() > epsn_x_max or self.bunch.epsn_y() > epsn_y_max
            ):
                stop = True
                print("Stop simulation due to emittance growth.")

        return stop

    def _finalize_multijob_mode(self):

        # save data for multijob operation and launch new job
        import h5py

        with h5py.File(
            "bunch_status_part%02d.h5" % (self.SimSt.present_simulation_part), "w"
        ) as fid:
            fid["bunch"] = self.piece_to_buffer(self.bunch)
        if not self.SimSt.first_run:
            os.system(
                "rm bunch_status_part%02d.h5"
                % (self.SimSt.present_simulation_part - 1)
            )
        self.SimSt.after_simulation()


def get_sim_instance(N_cores_pretend, id_pretend, init_sim_objects_auto=True):

    import PyPARIS.util as pu
    sim_instance = pu.get_sim_instance(Simulation(), N_cores_pretend, id_pretend,
        init_sim_objects_auto)

    return sim_instance


def get_serial_CPUring(init_sim_objects_auto=True):

    import PyPARIS.util as pu
    ring = pu.get_serial_CPUring(Simulation(), init_sim_objects_auto)

    return ring


class ParticleTrajectories(object):
    def __init__(self, n_record, n_turns):

        # prepare storage for particles coordinates
        self.x_i = np.empty((n_record, n_turns))
        self.xp_i = np.empty((n_record, n_turns))
        self.y_i = np.empty((n_record, n_turns))
        self.yp_i = np.empty((n_record, n_turns))
        self.z_i = np.empty((n_record, n_turns))
        self.i_turn = 0

    def dump(self, bunch):

        # id and momenta after track
        id_after = bunch.id
        x_after = bunch.x
        y_after = bunch.y
        z_after = bunch.z
        xp_after = bunch.xp
        yp_after = bunch.yp

        # sort id and momenta after track
        indsort = np.argsort(id_after)
        id_after = np.take(id_after, indsort)
        x_after = np.take(x_after, indsort)
        y_after = np.take(y_after, indsort)
        z_after = np.take(z_after, indsort)
        xp_after = np.take(xp_after, indsort)
        yp_after = np.take(yp_after, indsort)

        self.x_i[:, self.i_turn] = x_after
        self.xp_i[:, self.i_turn] = xp_after
        self.y_i[:, self.i_turn] = y_after
        self.yp_i[:, self.i_turn] = yp_after
        self.z_i[:, self.i_turn] = z_after

        self.i_turn += 1

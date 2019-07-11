### Some thoughts:
# - Input in the machine parameters
# - We introduce a cloud_dict
import scipy.io as sio
import os
import numpy as np
from numpy.random import rand
from scipy.constants import e as qe


class Ionization_Process(object):

    def __init__(self, pyecl_input_folder, process_name, process_definitions, cloud_dict):

        # Warn if target density doesn't correspond to density of gas ionization class?

        # Decide where to take into account the mass of the projectile (if not electrons, what do you need to do?)
        # - Use actual mass of projectile here (i.e. make sure that cross sections contain scaling to electron mass if needed)

        self.name = process_name
        print('Init process %s' % self.name)
        self.target_dens = process_definitions['target_density']
        self.E_eV_init = process_definitions['E_eV_init']

        #  Check that ionization product names correspond to existing clouds
        product_names = process_definitions['products']
        for product in product_names:
            assert product in cloud_dict.keys(), "Product name %s does not correspond to a defined cloud name."%(product)
        self.products = product_names

        # Read cross section file
        cross_section_file = process_definitions['cross_section']

        if os.path.isfile(pyecl_input_folder + '/' + cross_section_file):
            cross_section_file_path = pyecl_input_folder + '/' + cross_section_file
        elif os.path.isfile(pyecl_input_folder + '/' + cross_section_file + '.mat'):
            cross_section_file_path = pyecl_input_folder + '/' + cross_section_file + '.mat'
        else:
            cross_section_file_path = cross_section_file

        print('Cross-section from file %s' %cross_section_file_path)

        cross_section = sio.loadmat(cross_section_file_path)

        self.energy_eV = cross_section['energy_eV'].squeeze()
        self.sigma_cm2 = cross_section['cross_section_cm2'].squeeze()

        self.energy_eV_min = self.energy_eV.min()
        self.energy_eV_max = self.energy_eV.max()
        # Warn if minimum energy is not 0??

        # sey_diff is needed by the interp function
        # A 0 is appended because this last element is never needed but the array must have the correct shape
        self.sigma_cm2_diff = np.append(np.diff(self.sigma_cm2), 0.)

        flag_log = False

        # Check the energy step and define helpers for interp
        x_interp = self.energy_eV
        diff_x_interp = np.round(np.diff(x_interp), 3)
        delta_x_interp = diff_x_interp[0]
        x_interp_min = self.energy_eV_min

        if np.any(diff_x_interp != delta_x_interp):
            # Step not linear, check if logarithmic
            x_interp = np.log10(self.energy_eV)
            diff_x_interp = np.round(np.diff(x_interp), 3)
            delta_x_interp = diff_x_interp[0]
            x_interp_min = np.log10(self.energy_eV_min)

            if np.any(diff_x_interp != delta_x_interp):
                # Step neither linear nor logarithmic
                raise ValueError('Energy in cross section file must be equally spaced in linear or log scale.')
            else:
                flag_log = True

        self.delta_x_interp = delta_x_interp
        self.x_interp_min = x_interp_min
        self.flag_log = flag_log


    def get_sigma(self, energy_eV_proj):

        sigma_cm2_proj = energy_eV_proj * 0.

        # For now we set sigma = 0. both below and above energies in file...
        mask_below = (energy_eV_proj < self.energy_eV_min)
        mask_above = (energy_eV_proj > self.energy_eV_max)
        mask_interp = ~mask_below * ~mask_above

        if self.flag_log:
            x_interp_proj = np.log10(energy_eV_proj[mask_interp])
        else:
            x_interp_proj = energy_eV_proj[mask_interp]

        sigma_cm2_proj[mask_interp] = self._interp(x_interp_proj=x_interp_proj)

        # Return cross section in m2
        return sigma_cm2_proj * 1e-4 


    def _interp(self, x_interp_proj):
        """
        Linear interpolation of the energy - sigma curve.
        """
        index_float = (x_interp_proj - self.x_interp_min) / self.delta_x_interp
        index_remainder, index_int = np.modf(index_float)
        index_int = index_int.astype(int)

        return self.sigma_cm2[index_int] + index_remainder * self.sigma_cm2_diff[index_int]



class Cross_Ionization(object):

    def __init__(self, pyecl_input_folder, cross_ion_definitions, cloud_list):
        
        print('Initializing cross ionization.')

        # Make cloud dict from list
        self.cloud_dict = {}
        for cloud in cloud_list:
            self.cloud_dict.update({cloud.name : cloud})

        self.projectiles_dict = {}

        # Init projectiles
        for projectile in cross_ion_definitions.keys():
            print('Projectile %s:' %(projectile))
 
           # Check that projectile name corresponds to existing cloud
            assert projectile in self.cloud_dict.keys(), "Projectile name %s does not correspond to a defined cloud name."%(projectile)

            self.projectiles_dict.update({projectile : []})

            # Init processes
            for process_name in cross_ion_definitions[projectile].keys():
                process_definitions = cross_ion_definitions[projectile][process_name]
                process = Ionization_Process(pyecl_input_folder, process_name, process_definitions, self.cloud_dict)

                self.projectiles_dict[projectile].append(process)


    def generate(self, Dt):
        
        for projectile in self.projectiles_dict.keys():
            
            thiscloud_proj = self.cloud_dict[projectile]
            MP_e_proj = thiscloud_proj.MP_e
            N_proj = MP_e_proj.N_mp

            if N_proj > 0:

                # Calculate projectile energies
                vx_mp_proj = MP_e_proj.vx_mp[:N_proj]
                vy_mp_proj = MP_e_proj.vy_mp[:N_proj]
                vz_mp_proj = MP_e_proj.vz_mp[:N_proj]

                v_mp_proj = np.sqrt(vx_mp_proj * vx_mp_proj +
                                    vy_mp_proj * vy_mp_proj +
                                    vz_mp_proj * vz_mp_proj)

                E_eV_mp_proj = 0.5 * MP_e_proj.mass / qe * v_mp_proj * v_mp_proj


                for process in self.projectiles_dict[projectile]:

                    sigma_proj_MPs = process.get_sigma(energy_eV_proj=E_eV_mp_proj)

                    # Compute N_mp to add
                    DN_per_proj = sigma_proj_MPs * process.target_dens * v_mp_proj * Dt * MP_e_proj.nel_mp[:N_proj]

                    E_eV_init_gen = process.E_eV_init

                    for product in process.products:
                        thiscloud_gen = self.cloud_dict[product]
                        MP_e_gen = thiscloud_gen.MP_e
                        nel_mp_ref_gen = MP_e_gen.nel_mp_ref

                        # For now initialize generated MPs with velocity determined by input initial energy -
                        # similarly to gas ionization
                        v0_gen = -np.sqrt(2 * (E_eV_init_gen / 3.) * qe / MP_e_gen.mass)

                        ####################### Move to dedicated function #####################################
                        # N_mp_per_proj_float = DN_per_proj / MP_e_gen.nel_mp_ref

                        # N_mp_per_proj_int = np.floor(N_mp_per_proj_float)
                        # rest = N_mp_per_proj_float - N_mp_per_proj_int
                        # N_mp_per_proj_int = np.int_(N_mp_per_proj_int)

                        # N_mp_per_proj_int += np.int_(rand(N_proj) < rest)

                        # # For now initialize generated MPs with velocity determined by input initial energy -
                        # # similarly to gas ionization

                        # N_new_MPs = np.sum(N_mp_per_proj_int)
                        # print('Generating %d MPs in cloud %s' %(N_new_MPs, product))
                        # if N_new_MPs > 0:
                        #     # Compute MP_size
                        #     mask_gen = N_mp_per_proj_int > 0
                        #     N_mp_per_proj_int_masked = N_mp_per_proj_int[mask_gen]

                        #     nel_new_MPs_masked = np.zeros(np.sum(mask_gen))
                        #     nel_new_MPs_masked = DN_per_proj[mask_gen] / np.float_(N_mp_per_proj_int_masked)

                        #     nel_new_MPs = np.repeat(nel_new_MPs_masked, N_mp_per_proj_int_masked)

                        #     x_masked = MP_e_proj.x_mp[:N_proj][mask_gen]
                        #     y_masked = MP_e_proj.y_mp[:N_proj][mask_gen]
                        #     z_masked = MP_e_proj.z_mp[:N_proj][mask_gen]

                        #     x_new_MPs = np.repeat(x_masked, N_mp_per_proj_int_masked)
                        #     y_new_MPs = np.repeat(y_masked, N_mp_per_proj_int_masked)
                        #     z_new_MPs = np.repeat(z_masked, N_mp_per_proj_int_masked)

                        #     vx_new_MPs = np.zeros(N_new_MPs)
                        #     vy_new_MPs = np.zeros(N_new_MPs)
                        #     vz_new_MPs = np.zeros(N_new_MPs)

                        #     vx_new_MPs = v0_gen * (rand(N_new_MPs) - 0.5) # if you note a downwards polarization look here
                        #     vy_new_MPs = v0_gen * (rand(N_new_MPs) - 0.5)
                        #     vz_new_MPs = v0_gen * (rand(N_new_MPs) - 0.5)
                        ################### End move to dedicated function #####################################

                        (N_new_MPs, nel_new_MPs, x_new_MPs, y_new_MPs,
                         z_new_MPs, vx_new_MPs, vy_new_MPs,
                         vz_new_MPs) = get_new_mps(DN_per_proj=DN_per_proj,
                                                  N_mp_proj=N_proj,
                                                  nel_mp_ref_gen=nel_mp_ref_gen,
                                                  v0_gen=v0_gen)

                        if N_new_MPs > 0:
                            print('Generating %d MPs in cloud %s' %(N_new_MPs, product))
                            t_last_impact = -1
                            MP_e_gen.add_new_MPs(N_new_MPs, nel_new_MPs, x_new_MPs, y_new_MPs,
                                                 z_new_MPs, vx_new_MPs, vy_new_MPs, vz_new_MPs, t_last_impact)


    def get_new_mps(self, DN_per_proj, N_mp_proj, nel_mp_ref_gen, v0_gen):

        N_mp_per_proj_float = DN_per_proj / nel_mp_ref_gen

        N_mp_per_proj_int = np.floor(N_mp_per_proj_float)
        rest = N_mp_per_proj_float - N_mp_per_proj_int
        N_mp_per_proj_int = np.int_(N_mp_per_proj_int)
        N_mp_per_proj_int += np.int_(rand(N_mp_proj) < rest)

        N_new_MPs = np.sum(N_mp_per_proj_int)

        if N_new_MPs > 0:

            mask_gen = N_mp_per_proj_int > 0
            N_mp_per_proj_int_masked = N_mp_per_proj_int[mask_gen]

            nel_new_MPs_masked = np.zeros(np.sum(mask_gen))
            nel_new_MPs_masked = DN_per_proj[mask_gen] / np.float_(N_mp_per_proj_int_masked)

            nel_new_MPs = np.repeat(nel_new_MPs_masked, N_mp_per_proj_int_masked)

            x_masked = MP_e_proj.x_mp[:N_mp_proj][mask_gen]
            y_masked = MP_e_proj.y_mp[:N_mp_proj][mask_gen]
            z_masked = MP_e_proj.z_mp[:N_mp_proj][mask_gen]

            x_new_MPs = np.repeat(x_masked, N_mp_per_proj_int_masked)
            y_new_MPs = np.repeat(y_masked, N_mp_per_proj_int_masked)
            z_new_MPs = np.repeat(z_masked, N_mp_per_proj_int_masked)

            vx_new_MPs = np.zeros(N_new_MPs)
            vy_new_MPs = np.zeros(N_new_MPs)
            vz_new_MPs = np.zeros(N_new_MPs)

            vx_new_MPs = v0_gen * (rand(N_new_MPs) - 0.5)
            vy_new_MPs = v0_gen * (rand(N_new_MPs) - 0.5)
            vz_new_MPs = v0_gen * (rand(N_new_MPs) - 0.5)

        else:

            nel_new_MPs = np.array([])
            x_new_MPs = np.array([])
            y_new_MPs = np.array([])
            z_new_MPs = np.array([])
            vx_new_MPs = np.array([])
            vy_new_MPs = np.array([])
            vz_new_MPs = np.array([])

        return N_new_MPs, nel_new_MPs, x_new_MPs, y_new_MPs,\
               z_new_MPs, vx_new_MPs, vy_new_MPs, vz_new_MPs


    def extract_sigma(self, n_rep, energy_eV_test, Dt_test):

        sigmas = {}
        nel_mp = np.ones(n_rep)

        for projectile in self.projectiles_dict.keys():

            thiscloud_proj = self.cloud_dict[projectile]
            mass_proj = thiscloud_proj.MP_e.mass
            v_test = np.sqrt(2 * energy_eV_test * qe / mass_proj)
            # * np.ones_like(nel_mp)

            sigmas[projectile] = {}

            for process in self.projectiles_dict[projectile]:

                sigmas[projectile][process.name] = np.zeros(len(energy_eV_test))

                for i_ene, energy in enumerate(energy_eV_test):

                    energy_eV_array = energy * np.ones(nel_rep)
                    sigma_proj_MPs = process.get_sigma(energy_array)

                    # Compute N_mp to add
                    DN_per_proj = sigma_proj_MPs * process.target_dens * v_test * Dt_test * nel_mp

                    N_new_MPs, nel_new_MPs, x_new_MPs, y_new_MPs, z_new_MPs, vx_new_MPs, \
                        vy_new_MPs, vz_new_MPs = get_new_mps(DN_per_proj=DN_per_proj, N_mp_proj=N_proj,
                                                             nel_mp_ref_gen=MP_e_gen.nel_mp_ref, v0_gen=v0_gen)

                    sigmas[projectile][process.name][i_ene] = 

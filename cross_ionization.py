### Some thoughts:
# - Input in the machine parameters
# - We introduce a cloud_dict
import scipy.io as sio
import os
import numpy as np
from numpy.random import rand
from scipy.constants import e as qe

cross_ion_definitions = {
    'electrons': {
        'process e1': {'target_density': 1e10, 'products':['hydrogen', 'electrons'], 'cross_section':'hydrogen_from_electrons.mat', 'E_eV_init' : 0.1},
        'process e2': {'target_density': 1e10, 'products':['nitrogen', 'electrons'], 'cross_section':'nitrogen_from_electrons.mat', 'E_eV_init' : 0.1},
        },
    'nitrogen': {
        'process n1': {'target_density': 1e10, 'products':['hydrogen', 'electrons'], 'cross_section':'hydrogen_from_nitrogen.mat', 'E_eV_init' : 0.1},
        'process n2': {'target_density': 1e10, 'products':['nitrogen', 'electrons'], 'cross_section':'nitrogen_from_nitrogen.mat', 'E_eV_init' : 0.1},
        }
    }

class Ionization_Process(object):

    def __init__(self, pyecl_input_folder, process_name, process_definitions, cloud_dict):

        # Warn if target density doesn't correspond to density of gas ionization class?

        # Decide where to take into account the mass of the projectile (if not electrons, what do you need to do?)
        # - Use actual mass of projectile here (i.e. make sure that cross sections contain scaling to electron mass if needed)

        # Inspiration from sec_emission_model_from_file.py

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

        sigma_proj = energy_eV_proj * 0.

        mask_below = (energy_eV_proj < self.energy_eV_min)
        mask_above = (energy_eV_proj > self.energy_eV_max)
        mask_interp = ~mask_below * ~mask_above

        # For now we set sigma = 0. both below and above energies in file...

        if self.flag_log:
            x_interp_proj = np.log10(energy_eV_proj[mask_interp])
        else:
            x_interp_proj = energy_eV_proj[mask_interp]

        sigma_proj[mask_interp] = self._interp(x_interp_proj)

        return sigma_proj * 1e-4 # Go from cm2 to m2

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
        
        # Assert and warn as much as wanted:
        #  - Check that projectile names and ionization product names correspond to existing clouds

        # Inspiration from sec_emission_model_from_file.py

        print('Initializing cross ionization.')

        # Make cloud dict from list
        self.cloud_dict = {}
        for cloud in cloud_list:
            self.cloud_dict.update({cloud.name : cloud})

        self.projectiles_dict = {}

        for projectile in cross_ion_definitions.keys():
            print('Projectile %s:' %(projectile))
            assert projectile in self.cloud_dict.keys(), "Projectile name %s does not correspond to a defined cloud name."%(projectile)

            self.projectiles_dict.update({projectile : []})

            for process_name in cross_ion_definitions[projectile].keys():
                process_definitions = cross_ion_definitions[projectile][process_name]
                process = Ionization_Process(pyecl_input_folder, process_name, process_definitions, self.cloud_dict)

                self.projectiles_dict[projectile].append(process)

                ############### Moved to Process class ##################################

                # product_names = cross_ion_def_process['products']
                # for product in product_names:
                #     assert product in cloud_dict.keys(), "Product name %s does not correspond to a defined cloud name."%(product)

                # # Read cross section file
                # cross_section_file = cross_ion_def_process['cross_section']

                # if os.path.isfile(pyecl_input_folder + '/' + cross_section_file):
                #     cross_section_file_path = pyecl_input_folder + '/' + cc.filename_chm
                # elif os.path.isfile(pyecl_input_folder + '/' + cross_section_file + '.mat'):
                #     cross_section_file_path = pyecl_input_folder + '/' + cc.filename_chm + '.mat'
                # else:
                #     cross_section_file_path = cross_section_file

                # print('Cross-section from file %s' %cross_section_file_path)

                # cross_section = sio.loadmat(cross_section_file_path)

                # energy_eV = cross_section['energy_eV'].squeeze()
                # sigma_cm2 = cross_section['cross_section_cm2'].squeeze()

                # flag_log = False
                # diff_e = np.diff(energy_eV)
                # delta_e = diff_e[0]
                # if np.any(diff_e != delta_e):
                #     # Check if logarithmic
                #     energy_eV = np.log10(energy_eV)
                #     diff_e = np.diff(energy_eV)
                #     delta_e = diff_e[0]
                #     if np.any(diff_e != delta_e):
                #         raise ValueError('Energy in cross section file must be equally spaced in linear or log scale.')
                #     else:
                #         flag_log = True

                # cross_ion_def_process['energy_eV_sigma'] = energy_eV
                # cross_ion_def_process['sigma_cm2'] = sigma_cm2
                # cross_ion_def_process['flag_log'] = flag_log

                ############### End moved to Process class ##################################


    def generate(self, Dt):
        
        for projectile in self.projectiles_dict.keys():
            
            thiscloud_proj = self.cloud_dict[projectile]
            MP_e_proj = thiscloud_proj.MP_e
            N_proj = MP_e_proj.N_mp
 
            if N_proj > 0:

                vx_mp_proj = MP_e_proj.vx_mp[:N_proj]
                vy_mp_proj = MP_e_proj.vy_mp[:N_proj]
                vz_mp_proj = MP_e_proj.vz_mp[:N_proj]

                v_mp_proj = np.sqrt(vx_mp_proj * vx_mp_proj +
                                    vy_mp_proj * vy_mp_proj +
                                    vz_mp_proj * vz_mp_proj)

                E_eV_mp_proj = 0.5 * MP_e_proj.mass / qe * v_mp_proj * v_mp_proj

                # Evaluate from the speeds in the corresponding clouds
                """
                Something like:
                v_impact_mod = np.sqrt(vx_impact * vx_impact + vy_impact * vy_impact + vz_impact * vz_impact)
                E_impact_eV = 0.5 * MP_e.mass / qe * v_impact_mod * v_impact_mod
                """

                for process in self.projectiles_dict[projectile]:

                    # cross_ion_def_process = cross_ion_def_projectile[process]

                    # dens = cross_ion_def_process['target_density']
                    # sigma_cm2 = cross_ion_def_process['sigma_cm2']
                    # energy_eV_sigma = cross_ion_def_process['energy_eV_sigma']
                    # product_list = cross_ion_def_process['products']

                    # Interpolate cross section
                    # Inspiration from sec_emission_model_from_file.py
                    
                    sigma_proj_MPs = process.get_sigma(E_eV_mp_proj)

                    # Compute N_mp to add
                    DN_per_proj = sigma_proj_MPs * process.target_dens * v_mp_proj * Dt * MP_e_proj.nel_mp[:N_proj]

                    E_eV_init_gen = process.E_eV_init

                    for product in process.products:
                        thiscloud_gen = self.cloud_dict[product]
                        MP_e_gen = thiscloud_gen.MP_e

                        N_mp_per_proj_float = DN_per_proj / MP_e_gen.nel_mp_ref

                        N_mp_per_proj_int = np.floor(N_mp_per_proj_float)
                        rest = N_mp_per_proj_float - N_mp_per_proj_int
                        N_mp_per_proj_int = np.int_(N_mp_per_proj_int)

                        N_mp_per_proj_int += np.int_(rand(N_proj) < rest)

                        # For now initialize generated MPs with velocity determined by input initial energy -
                        # similarly to gas ionization
                        v0 = -np.sqrt(2 * (E_eV_init_gen / 3.) * qe / MP_e_gen.mass)

                        N_new_MPs = np.sum(N_mp_per_proj_int)
                        print('Generating %d MPs in cloud %s' %(N_new_MPs, product))
                        if N_new_MPs > 0:
                            # Compute MP_size
                            mask_gen = N_mp_per_proj_int > 0
                            N_mp_per_proj_int_masked = N_mp_per_proj_int[mask_gen]

                            nel_new_MPs_masked = np.zeros(np.sum(mask_gen))
                            nel_new_MPs_masked = DN_per_proj[mask_gen] / np.float_(N_mp_per_proj_int_masked)

                            nel_new_MPs = np.repeat(nel_new_MPs_masked, N_mp_per_proj_int_masked)

                            x_masked = MP_e_proj.x_mp[:N_proj][mask_gen]
                            y_masked = MP_e_proj.y_mp[:N_proj][mask_gen]
                            z_masked = MP_e_proj.z_mp[:N_proj][mask_gen]

                            x_new_MPs = np.repeat(x_masked, N_mp_per_proj_int_masked)
                            y_new_MPs = np.repeat(y_masked, N_mp_per_proj_int_masked)
                            z_new_MPs = np.repeat(z_masked, N_mp_per_proj_int_masked)

                            vx_new_MPs = np.zeros(N_new_MPs)
                            vy_new_MPs = np.zeros(N_new_MPs)
                            vz_new_MPs = np.zeros(N_new_MPs)

                            vx_new_MPs = v0 * (rand(N_new_MPs) - 0.5) # if you note a towards down polarization look here
                            vy_new_MPs = v0 * (rand(N_new_MPs) - 0.5)
                            vz_new_MPs = v0 * (rand(N_new_MPs) - 0.5)

                            t_last_impact = -1
                            MP_e_gen.add_new_MPs(N_new_MPs, nel_new_MPs, x_new_MPs, y_new_MPs,
                                                 z_new_MPs, vx_new_MPs, vy_new_MPs, vz_new_MPs, t_last_impact)






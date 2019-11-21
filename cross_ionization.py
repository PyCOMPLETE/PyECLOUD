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
#                           Lorenzo Giacomel
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

import scipy.io as sio
import os
import numpy as np
from numpy.random import rand
from scipy.constants import e as qe


class Ionization_Process(object):

    def __init__(self, pyecl_input_folder, process_name, process_definitions, cloud_dict, target_area):

        # Warn if target density doesn't correspond to density of gas ionization class?

        self.name = process_name
        print('Init process %s' % self.name)

        self.target_dens = process_definitions['target_density']
        print('Target density = %.2e' %(self.target_dens))
        self.last_reported_target_dens = self.target_dens

        self.target_area = target_area
        self.N_target = self.target_dens * target_area

        self.E_eV_init = process_definitions['E_eV_init']

        if 'extract_sigma' in process_definitions.keys():
            self.extract_sigma = process_definitions['extract_sigma']
        else:
            self.extract_sigma = True

        if 'generate_equally' in process_definitions.keys():
            self.generate_equally = process_definitions['generate_equally']
        else:
            self.generate_equally = False

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

        if self.extract_sigma:
            self.extract_sigma_path = cross_section_file_path.split('.mat')[0]
            self.extract_sigma_path += '_extracted.mat'
        else:
            self.extract_sigma_path = None

        self.energy_eV = cross_section['energy_eV'].squeeze()
        self.sigma_cm2 = cross_section['cross_section_cm2'].squeeze()

        # Check the energy step and define helpers for interp
        self.energy_eV_min = self.energy_eV.min()
        self.energy_eV_max = self.energy_eV.max()

        self.sigma_cm2_diff = np.append(np.diff(self.sigma_cm2), 0.)
        # A 0 is appended to give the array the correct shape

        flag_log = False

        ndec_round_x = 8
        x_interp = self.energy_eV
        diff_x_interp = np.round(np.diff(x_interp), ndec_round_x)
        delta_x_interp = diff_x_interp[0]
        x_interp_min = self.energy_eV_min

        if np.any(diff_x_interp != delta_x_interp):
            # Step not linear, check if logarithmic
            x_interp = np.log10(self.energy_eV)
            diff_x_interp = np.round(np.diff(x_interp), ndec_round_x)
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


    def generate(self, Dt, cloud_dict, mass_proj, N_proj, nel_mp_proj,
                 x_proj, y_proj, z_proj, v_mp_proj, flag_generate=True):

        E_eV_mp_proj = 0.5 * mass_proj / qe * v_mp_proj * v_mp_proj

        # Get sigma
        sigma_mp_proj = self.get_sigma(energy_eV_proj=E_eV_mp_proj)

        DN_per_proj = sigma_mp_proj * self.target_dens * v_mp_proj * Dt * nel_mp_proj

        N_proj = len(nel_mp_proj)

        # Calculate remaining density
        if flag_generate:
            DN_target = np.sum(DN_per_proj)
            self.N_target =  np.round(self.N_target - DN_target, 3)
            self.target_dens = self.N_target / self.target_area

            if self.target_dens < 0.1 * self.last_reported_target_dens:
                print('Cross-ionization process %s target density = %.2e' %(self.name, self.target_dens))
                self.last_reported_target_dens = self.target_dens

        new_mp_info = {}

        if self.generate_equally:
            # Calculate average product nel_mp_ref
            nel_mp_ref_products = 0.
            N_products = len(self.products)

            for product in self.products:
                thiscloud_gen = cloud_dict[product]
                MP_e_gen = thiscloud_gen.MP_e
                nel_mp_ref_products += MP_e_gen.nel_mp_ref / N_products

            # Compute N_mp to add (the same for all products)
            N_mp_per_proj_float = DN_per_proj / nel_mp_ref_products
            N_mp_per_proj_int = np.floor(N_mp_per_proj_float)
            rest = N_mp_per_proj_float - N_mp_per_proj_int
            N_mp_per_proj_int = np.atleast_1d(np.int_(N_mp_per_proj_int))
            N_mp_per_proj_int += np.atleast_1d(np.int_(rand(N_proj) < rest))

            N_new_MPs = np.sum(N_mp_per_proj_int)

        for product in self.products:

            new_mp_info[product] = {}

            thiscloud_gen = cloud_dict[product]
            MP_e_gen = thiscloud_gen.MP_e
            mass_gen = MP_e_gen.mass

            # Initialize generated MPs with energy defined by user 
            v0_gen = np.sqrt(2 * (self.E_eV_init / 3.) * qe / mass_gen)

            if self.generate_equally:
                nel_mp_ref_gen = nel_mp_ref_products
            else:
                nel_mp_ref_gen = MP_e_gen.nel_mp_ref

                # Compute N_mp to add (different for each product)
                N_mp_per_proj_float = DN_per_proj / nel_mp_ref_gen
                N_mp_per_proj_int = np.floor(N_mp_per_proj_float)
                rest = N_mp_per_proj_float - N_mp_per_proj_int
                N_mp_per_proj_int = np.atleast_1d(np.int_(N_mp_per_proj_int))
                N_mp_per_proj_int += np.atleast_1d(np.int_(rand(N_proj) < rest))

                N_new_MPs = np.sum(N_mp_per_proj_int)

            if N_new_MPs > 0:
                mask_gen = N_mp_per_proj_int > 0
                N_mp_per_proj_int_masked = N_mp_per_proj_int[mask_gen]

                nel_new_MPs_masked = np.ones(np.sum(mask_gen)) * nel_mp_ref_gen

                nel_new_MPs = np.repeat(nel_new_MPs_masked, N_mp_per_proj_int_masked)

                x_masked = x_proj[mask_gen]
                y_masked = y_proj[mask_gen]
                z_masked = z_proj[mask_gen]

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

            new_mp_info[product]['N_new_MPs'] = N_new_MPs
            new_mp_info[product]['nel_new_MPs'] = nel_new_MPs
            new_mp_info[product]['x_new_MPs'] = x_new_MPs
            new_mp_info[product]['y_new_MPs'] = y_new_MPs
            new_mp_info[product]['z_new_MPs'] = z_new_MPs
            new_mp_info[product]['vx_new_MPs'] = vx_new_MPs
            new_mp_info[product]['vy_new_MPs'] = vy_new_MPs
            new_mp_info[product]['vz_new_MPs'] = vz_new_MPs

        return new_mp_info, np.sum(DN_per_proj)


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

    def __init__(self, pyecl_input_folder, cross_ion_definitions, cloud_list,
                 chamber_area, n_rep_test=10000, Dt_test=25e-11,
                 energy_eV_test=np.logspace(np.log10(1.), np.log10(25000.), num=5000)):
        
        print('Initializing cross ionization.')

        # Make cloud dict from list
        cloud_dict = {}
        for cloud in cloud_list:
            cloud_dict.update({cloud.name : cloud})

        self.projectiles_dict = {}
        self.products = []

        # Init projectiles and make list of products
        for projectile in cross_ion_definitions.keys():
            print('Projectile %s:' %(projectile))
 
           # Check that projectile name corresponds to existing cloud
            assert projectile in cloud_dict.keys(), "Projectile name %s does not correspond to a defined cloud name."%(projectile)

            self.projectiles_dict.update({projectile : []})

            # Init processes
            for process_name in cross_ion_definitions[projectile].keys():
                process_definitions = cross_ion_definitions[projectile][process_name]
                process = Ionization_Process(pyecl_input_folder, process_name,
                                             process_definitions, cloud_dict, chamber_area)

                self.projectiles_dict[projectile].append(process)

                for product in process.products:
                    if product not in self.products:
                        self.products.append(product)

        # Extract sigma curves for consistency checks
        self._extract_sigma(Dt=Dt_test, cloud_dict=cloud_dict,
                            n_rep=n_rep_test, energy_eV=energy_eV_test)

        # Initialize dictionary for quantities to save
        self.nel_cross_ion = {}
        self.N_mp_cross_ion = {}
        self.DN_proj = {}
        for cloud in cloud_list:
            self.nel_cross_ion[cloud.name] = 0.
            self.N_mp_cross_ion[cloud.name] = 0
            self.DN_proj[cloud.name] = 0.


    def generate(self, Dt, cloud_list):
        
        # Make cloud dict from list
        cloud_dict = {}
        for cloud in cloud_list:
            cloud_dict.update({cloud.name : cloud})

        new_mps_to_gen = self._init_new_mp_dict(self.products)

        for projectile in self.projectiles_dict.keys():
            thiscloud = cloud_dict[projectile]
            MP_e = thiscloud.MP_e
            N_mp = MP_e.N_mp
            mass = MP_e.mass

            if N_mp > 0:

                nel_mp = MP_e.nel_mp[:N_mp]

                x_mp = MP_e.x_mp[:N_mp]
                y_mp = MP_e.y_mp[:N_mp]
                z_mp = MP_e.z_mp[:N_mp]

                vx_mp = MP_e.vx_mp[:N_mp]
                vy_mp = MP_e.vy_mp[:N_mp]
                vz_mp = MP_e.vz_mp[:N_mp]

                v_mp = np.sqrt(vx_mp * vx_mp +
                               vy_mp * vy_mp +
                               vz_mp * vz_mp)

                for process in self.projectiles_dict[projectile]:

                    mp_info_from_proc, DN_proj = process.generate(Dt=Dt,
                                                         cloud_dict=cloud_dict,
                                                         mass_proj=mass,
                                                         N_proj=N_mp,
                                                         nel_mp_proj=nel_mp,
                                                         x_proj=x_mp,
                                                         y_proj=y_mp,
                                                         z_proj=z_mp,
                                                         v_mp_proj=v_mp)

                    for product in process.products:
                        self._add_to_mp_dict(new_mps_to_gen[product],
                                             mp_info_from_proc[product])
                        self.DN_proj[product] += DN_proj

        t_last_impact = -1
        for thiscloud in cloud_list:
            if thiscloud.name in self.products:
                MP_e = thiscloud.MP_e
                new_mps = new_mps_to_gen[thiscloud.name]

                if new_mps['N_new_MPs'] > 0:
                    MP_e.add_new_MPs(new_mps['N_new_MPs'],
                                     new_mps['nel_new_MPs'],
                                     new_mps['x_new_MPs'],
                                     new_mps['y_new_MPs'],
                                     new_mps['z_new_MPs'],
                                     new_mps['vx_new_MPs'],
                                     new_mps['vy_new_MPs'],
                                     new_mps['vz_new_MPs'],
                                     t_last_impact)

                # Add to saved data
                self.nel_cross_ion[thiscloud.name] += np.sum(new_mps['nel_new_MPs'])
                self.N_mp_cross_ion[thiscloud.name] += new_mps['N_new_MPs']
            else:
                self.nel_cross_ion[thiscloud.name] += 0.
                self.N_mp_cross_ion[thiscloud.name] += 0


    def save_cross_ion_data(self, cloud_name):

        thiscloud_nel_cross_ion = self.nel_cross_ion[cloud_name]
        thiscloud_N_mp_cross_ion = self.N_mp_cross_ion[cloud_name]
        thiscloud_DN_proj = self.DN_proj[cloud_name]
        self.nel_cross_ion[cloud_name] = 0.
        self.N_mp_cross_ion[cloud_name] = 0.
        self.DN_proj[cloud_name] = 0.

        return thiscloud_nel_cross_ion, thiscloud_N_mp_cross_ion, thiscloud_DN_proj


    def _extract_sigma(self, Dt, cloud_dict, n_rep, energy_eV):

        v0 = 0.
        N_ene = len(energy_eV)

        N_mp = n_rep

        x_mp = np.zeros(n_rep)
        y_mp = np.zeros(n_rep)
        z_mp = np.zeros(n_rep)

        for projectile in self.projectiles_dict.keys():

            thiscloud = cloud_dict[projectile]
            mass = thiscloud.MP_e.mass
            nel_mp = np.ones(n_rep) * thiscloud.MP_e.nel_mp_ref

            v_test = np.sqrt(2 * energy_eV * qe / mass)

            for process in self.projectiles_dict[projectile]:

                if process.extract_sigma:

                    print('Extracting cross section for process %s' %process.name )

                    save_dict = {}
                    save_dict['energy_eV'] = energy_eV
                    save_dict['sigma_cm2_interp'] = np.zeros(len(energy_eV))
                    for product in process.products:
                        this_sigma_name = 'sigma_cm2_sampled_%s' %(product)
                        save_dict[this_sigma_name] = np.zeros(len(energy_eV))

                    for i_ene, energy in enumerate(energy_eV):

                        if np.mod(i_ene, N_ene / 10) == 0:
                            print ('Extracting sigma %.0f'%(float(i_ene) / float(N_ene) * 100) + """%""")

                        # Test process.get_sigma()
                        sigma_m2 = process.get_sigma(np.array([energy]))
                        save_dict['sigma_cm2_interp'][i_ene] = sigma_m2 * 1e4

                        # Test process.generate()
                        v_ene = v_test[i_ene]
                        v_mp = v_ene * np.ones(n_rep)

                        mp_info_from_proc, _ = process.generate(Dt, cloud_dict=cloud_dict,
                                                             mass_proj=mass, N_proj=N_mp,
                                                             nel_mp_proj=nel_mp, x_proj=x_mp,
                                                             y_proj=y_mp, z_proj=z_mp,
                                                             v_mp_proj=v_mp, flag_generate=False)

                        for product in mp_info_from_proc.keys():
                            DN_gen = np.sum(mp_info_from_proc[product]['nel_new_MPs'])
                            if v_ene > 0:
                                sigma_m2_est = DN_gen / process.target_dens / v_ene / Dt / np.sum(nel_mp)
                            else:
                                sigma_m2_est = 0.
                            this_sigma_name = 'sigma_cm2_sampled_%s' %(product)
                            save_dict[this_sigma_name][i_ene] = sigma_m2_est * 1e4

                    sio.savemat(process.extract_sigma_path, save_dict, oned_as='row')
                    print('Saved extracted cross section as %s' %process.extract_sigma_path)


    def _init_new_mp_dict(self, products):
        # Init new MP dictionary for products
        new_mp_dict = {}
        mp_dict_keys = ['N_new_MPs', 'nel_new_MPs', 'x_new_MPs', 'y_new_MPs',
                        'z_new_MPs', 'vx_new_MPs', 'vy_new_MPs', 'vz_new_MPs']
        for product in products:
            new_mp_dict[product] = {}
            for key in mp_dict_keys:
                if key == 'N_new_MPs':
                    new_mp_dict[product][key] = 0
                else:
                    new_mp_dict[product][key] = np.array([])

        return new_mp_dict


    def _add_to_mp_dict(self, mp_dict, dict_to_add):
        #sum_dict = {}
        for key in mp_dict.keys():
            if key == 'N_new_MPs':
                mp_dict[key] += dict_to_add['N_new_MPs']
            else:
                mp_dict[key] = np.append(mp_dict[key], dict_to_add[key])


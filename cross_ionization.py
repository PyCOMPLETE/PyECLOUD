### Some thoughts:
# - Input in the machine parameters
# - We introduce a cloud_dict
import scipy.io as sio
import os
from numpy.random import rand

cross_ion_definitions = {
    'electrons': {
        'process e1': {'target_density': 1e10, 'products':['hydrogen', 'electrons'], 'cross_section':'hydrogen_from_electrons.mat'},
        'process e2': {'target_density': 1e10, 'products':['nitrogen', 'electrons'], 'cross_section':'nitrogen_from_electrons.mat'},
        } 
    'nitrogen': {
        'process n1': {'target_density': 1e10, 'products':['hydrogen', 'electrons'], 'cross_section':'hydrogen_from_nitrogen.mat'},
        'process n2': {'target_density': 1e10, 'products':['nitrogen', 'electrons'], 'cross_section':'nitrogen_from_nitrogen.mat'},
        }
    }

 

class Cross_Ionization(object):

    def __init__(self, cross_ion_definitions, cloud_dict):
        
        self.cross_ion_definitions = cross_ion_definitions

        # Read the files here
        # Assert and warn as much as wanted:
        #  - Check that projectile names and ionization product names correspond to existing clouds
        # Update dictionary with data extracted from file (logE, sigma)

        # Decide where to take into account the mass of the projectile (if not electrons, what do you need to do?)
        # - Use actual mass of projectile here (i.e. make sure that cross sections contain scaling to electron mass if needed)

        # Inspiration from sec_emission_model_from_file.py

        print('Initializing cross ionization.')

        for projectile in self.cross_ion_definitions.keys():

            print('Projectile %s' %(projectile))

            assert projectile in cloud_dict.keys(), "Projectile name %s does not correspond to a defined cloud name."%(projectile)

            cross_ion_def_projectile = self.cross_ion_definitions[projectile]

            for process in cross_ion_def_projectile.keys():

                # Warn if target density doesn't correspond to density of gas ionization class?
                # dens = cross_ion_def_projectile[process]['target_density']

                product_names = cross_ion_def_projectile[process]['products']
                for product in product_names:
                    assert product in cloud_dict.keys(), "Product name %s does not correspond to a defined cloud name."%(product)

                # Read cross section file
                cross_section_file = cross_ion_def_projectile[process]['cross_section']

                if os.path.isfile(pyecl_input_folder + '/' + cross_section_file):
                    cross_section_file_path = pyecl_input_folder + '/' + cc.filename_chm
                elif os.path.isfile(pyecl_input_folder + '/' + cross_section_file + '.mat'):
                    cross_section_file_path = pyecl_input_folder + '/' + cc.filename_chm + '.mat'
                else:
                    cross_section_file_path = cross_section_file

                print('Cross-section from file %s' %cross_section_file_path)

                cross_section = sio.loadmat(cross_section_file_path)

                energy_eV_sigma = cross_section['energy_eV'].squeeze()
                sigma_cm2 = cross_section['cross_section_cm2'].squeeze()

                cross_ion_def_projectile[process]['energy_eV_sigma'] = energy_eV_sigma
                cross_ion_def_projectile[process]['sigma_cm2'] = sigma_cm2

    def generate(self, cloud_dict, Dt):
        
        for projectile in self.cross_ion_definitions.keys():
            
            this_cloud_proj = cloud_dict[projectile]
            MP_e_proj = this_cloud_proj.MP_e
            N_proj = MP_e_proj.N_mp
 
            if N_proj > 0:

                vx_mp_proj = MP_e_proj.vx_mp[0:N_proj]
                vy_mp_proj = MP_e_proj.vy_mp[0:N_proj]
                vz_mp_proj = MP_e_proj.vz_mp[0:N_proj]

                v_mp_proj = np.sqrt(vx_mp_proj * vx_mp_proj +
                                    vy_mp_proj * vy_mp_proj +
                                    vz_mp_proj * vz_mp_proj)

                E_mp_proj_eV = 0.5 * MP_e_proj.mass / qe * v_mp_proj * v_mp_proj

                # Evaluate from the speeds in the corresponding clouds
                """
                Something like:
                v_impact_mod = np.sqrt(vx_impact * vx_impact + vy_impact * vy_impact + vz_impact * vz_impact)
                E_impact_eV = 0.5 * MP_e.mass / qe * v_impact_mod * v_impact_mod
                """

                cross_ion_def_projectile = self.cross_ion_definitions[projectile]

                for process in cross_ion_def_projectile.keys():

                    dens = cross_ion_def_projectile[process]['target_density']
                    sigma_cm2 = cross_ion_def_projectile[process]['sigma_cm2']
                    logE_sigma = cross_ion_def_projectile[process]['energy_eV_sigma']
                    product_list = cross_ion_def_projectile[process]['products']

                    # Interpolate cross section
                    # Inspiration from sec_emission_model_from_file.py
                    
                    sigma_proj_MPs = 
                    
                    for cld_gen_name in generate_list:
                        this_cloud_gen =  cloud_dict[cld_gen_name]
                        MP_e_gen = this_cloud_gen.MP_e
                        
                        # Compute N_mp to add
                        DN_per_proj = sigma_proj_MPs*dens*v_proj_MPs*Dt*MP_e_proj.nel_mp[:N_proj]

                        N_mp_per_proj_float = DN_per_proj/MP_e_gen.nel_mp_ref

                        N_mp_per_proj_int = np.floor(N_mp_per_proj_float)
                        rest = N_mp_per_proj_float - N_mp_per_proj_int
                        N_mp_per_proj_int = np.int_(N_mp_per_proj_int)

                        N_mp_per_proj_int += np.int_(rand(N_proj) < rest)

                        if np.sum(N_mp_per_proj_int)>0:
                            # Compute MP_size
                            mask_gen = N_mp_per_proj_int>0
                            N_mp_per_proj_int_masked = N_mp_per_proj_int[mask_gen]

                            nel_mp_new =  np.zeros(np.sum(mask_gen))
                            nel_mp_new = DN_per_proj[mask_gen]/np.float_(N_mp_per_proj_int_masked)

                            x_masked = MP_e_proj...
                            y_masked = 
                            z_masked = 


                            x_new_MPs = np.repeat(x_masked, N_mp_per_proj_int_masked)
                            y_new_MPs = np.repeat(y_masked, N_mp_per_proj_int_masked)
                            z_new_MPs = np.repeat(z_masked, N_mp_per_proj_int_masked)

                            MP_e.add_new_MPs(N_new_MPs, nel_new_MPs, x_new_MPs, y_new_MPs,
                                             z_new_MPs,vx_new_MPs, vy_new_MPs, vz_new_MPs)










### Some thoughts:
# - Input in the machine parameters
# - We introduce a cloud_dict
from numpy.random import rand

cross_ion_definitions = {
    'electrons': {
        'process e1': {'target_density': 1e10, 'generate':['hydrogen', 'electrons'], 'cross_section':'hydrogen_from_electrons.mat'},
        'process e2': {'target_density': 1e10, 'generate':['nitrogen', 'electrons'], 'cross_section':'nitrogen_from_electrons.mat'},
        } 
    'nitrogen': {
        'process n1': {'target_density': 1e10, 'generate':['hydrogen', 'electrons'], 'cross_section':'hydrogen_from_nitrogen.mat'},
        'process n2': {'target_density': 1e10, 'generate':['nitrogen', 'electrons'], 'cross_section':'nitrogen_from_nitrogen.mat'},
        }
    }

 

class Cross_Ionization(object):

    def __init__(self, cross_ion_definitions, cloud_dict):
        
        self.cross_ion_definitions = cross_ion_definitions
        
        # Read the files here
        # Assert and warn as much as wanted
        # Update dictionary with data extracted from file (logE, sigma)

        # Inspiration from sec_emission_model_from_file.py


    def generate(self, cloud_dict, Dt):
        
        for projectile in cross_ion_definitions.keys():
            
            this_cloud_proj = cloud_dict[projectile]
            MP_e_proj = this_cloud_proj.MP_e
            N_proj = MP_e_proj.N_mp
 
            
            v_proj_MPs = # Evaluate from the speeds in the corresponding clouds
            """
            Something like:
            v_impact_mod = np.sqrt(vx_impact * vx_impact + vy_impact * vy_impact + vz_impact * vz_impact)
            E_impact_eV = 0.5 * MP_e.mass / qe * v_impact_mod * v_impact_mod
            """

            cross_ion_def_projectile = self.cross_ion_definitions[projectile]

            for process in cross_ion_def_projectile.keys():
                dens = cross_ion_def_projectile[process]['target_density']
                sigma = cross_ion_def_projectile[process]['sigma']
                logE_sigma = cross_ion_def_projectile[process]['logE_sigma']
                generate_list = cross_ion_def_projectile[process]['generate']

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

                                            MP_e.add_new_MPs(N_new_MPs, nel_new_MPs, x_new_MPs, y_new_MPs, z_new_MPs,
                                     vx_new_MPs, vy_new_MPs, vz_new_MPs)                       
                    









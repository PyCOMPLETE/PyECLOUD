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
#                   PyECLOUD Version 4.30                     
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
from scipy.constants import c
import time
from PyEC4PyHT import Ecloud
from gas_ionization_class import residual_gas_ionization


class MP_light(object):
    pass



class Ecloud_fastion(Ecloud):

    #@profile	
    def track(self, beam):

        start_time = time.mktime(time.localtime())
        self._reinitialize()

        MP_e = self.MP_e
        dynamics = self.dynamics
        impact_man = self.impact_man
        spacech_ele = self.spacech_ele

        if hasattr(beam.particlenumber_per_mp, '__iter__'):
            raise ValueError('ecloud module assumes same size for all beam MPs')
        
        if self.flag_clean_slices:
            beam.clean_slices()
        
        slices = beam.get_slices(self.slicer)
        self.slicer.add_statistics(sliceset=slices, beam=beam, statistics=True)

        # Only track over slices with particles
        filled_slices = np.where(slices.n_macroparticles_per_slice > 0)[0]
        
        for i in filled_slices[::-1]:

            # select particles in the bunch
            ix = slices.particle_indices_of_slice(i)
            
            # slice size and time step
            dz = (slices.z_bins[i + 1] - slices.z_bins[i]) # in this case it is the bucket length
            dt = dz / (beam.beta * c)
            
            # define substep
            if dt > self.Dt_ref:
                N_sub_steps = int(np.round(dt / self.Dt_ref))
            else:
                N_sub_steps = 1

            Dt_substep = dt / N_sub_steps
            # print Dt_substep, N_sub_steps, dt
            N_mp_bunch = slices.n_macroparticles_per_slice[i]
            
            if N_mp_bunch > 0:
                # beam field 
                # print 'Bunch', i
                MP_p = MP_light()
                MP_p.x_mp = beam.x[ix]
                MP_p.y_mp = beam.y[ix]
                MP_p.nel_mp = beam.x[ix] * 0. + beam.particlenumber_per_mp
                MP_p.N_mp = N_mp_bunch
                MP_p.charge = beam.charge


                if self.gas_ion_flag == 1:
                    Np_bunch = N_mp_bunch * slices.particlenumber_per_mp
                    dz_bunch = slices.slice_widths[i]
                    lambda_bunch = Np_bunch
                    dt_bunch = 1 / c
                    MP_e = self.gas_ionization.generate(MP_e=MP_e, lambda_t=lambda_bunch, Dt=dt_bunch, sigmax=slices.sigma_x[i], 
                                                    sigmay=slices.sigma_y[i], x_beam_pos=slices.mean_x[i], y_beam_pos=slices.mean_y[i])


                # compute beam field
                spacech_ele.recompute_spchg_efield(MP_p)

                # gather beam field to electrons
                Ex_n_beam, Ey_n_beam = spacech_ele.get_sc_eletric_field(MP_e) # since we do not divide the charge by the dz, this is already integrated over the bucket length

                # compute cloud field
                spacech_ele.recompute_spchg_efield(MP_e) 
                
                # gather cloud field to beam
                Ex_sc_p, Ey_sc_p = spacech_ele.get_sc_eletric_field(MP_p) 

                # kick cloud particles
                MP_e.vx_mp[:MP_e.N_mp] += Ex_n_beam * MP_e.charge / MP_e.mass / c
                MP_e.vy_mp[:MP_e.N_mp] += Ey_n_beam * MP_e.charge / MP_e.mass / c
                
                # kick beam particles
                fact_kick = beam.charge / (beam.mass * beam.beta * beam.beta * beam.gamma * c * c) * self.L_ecloud
                beam.xp[ix] += fact_kick * Ex_sc_p
                beam.yp[ix] += fact_kick * Ey_sc_p
            

            # Total electric field on electrons
            Ex_n = MP_e.vx_mp[:MP_e.N_mp] * 0.
            Ey_n = Ex_n

            # save position before motion step
            old_pos = MP_e.get_positions()
            
            # motion electrons
            MP_e = dynamics.stepcustomDt(MP_e, Ex_n,Ey_n, Dt_substep=Dt_substep, N_sub_steps=N_sub_steps)
            
            # impacts: backtracking and secondary emission
            MP_e = impact_man.backtrack_and_second_emiss(old_pos, MP_e)
            
                        
            if self.save_ele_distributions_last_track:
                self.rho_ele_last_track.append(spacech_ele.rho.copy())
                #print 'Here'
            
            if self.save_ele_potential_and_field:
                self.phi_ele_last_track.append(spacech_ele.phi.copy())
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
                

        if self.beam_monitor != None:
            self.beam_monitor.dump(beam)

        self._finalize()

        stop_time = time.mktime(time.localtime())
        print 'Done track in ', (stop_time-start_time), 's'
        #print self.N_MP_last_track
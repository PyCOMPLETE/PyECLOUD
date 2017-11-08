import numpy as np
from scipy.constants import c, e as qe

from . import input_parameters_format_specification as inp_spec
from . import myloadmat_to_obj as mlm

class beam_descr_from_fil:
    def __init__(self, beamfilename, betafx_from_mach_parms_file, Dx_from_mach_parms_file,
                 betafy_from_mach_parms_file, Dy_from_mach_parms_file):

        print('Parsing beam file: %s' % beamfilename)

        # Parse beam input file
        beam_beam = inp_spec.import_module_from_file('beam_beam', beamfilename)
        
        # Verify validity of provided module w.r.t. parameters_dict
        inp_spec.assert_module_has_parameters(beam_beam, 'beam_beam')
        
        # Create config_dict with all allowed beam parameters (not specified are set to default)
        config_dict = {}
        inp_spec.update_config_dict(config_dict, beam_beam, 'beam_beam')

        cc = mlm.obj_from_dict(config_dict)


        # keep beackwards compatibility (optics in the machine params file)
        if cc.Dx is None:
            cc.Dx = Dx_from_mach_parms_file
        if cc.Dy is None:
            cc.Dy = Dy_from_mach_parms_file
        if cc.betafx is None:
            cc.betafx = betafx_from_mach_parms_file
        if cc.betafy is None:
            cc.betafy = betafy_from_mach_parms_file

        energy_J = cc.energy_eV * qe
        gamma_rel= energy_J/(cc.m0_part*c**2)
        beta_rel = np.sqrt(1-1/(gamma_rel*gamma_rel))
        beta_rel = beta_rel.real

        if (cc.sigmax==-1) and ((cc.betafx is None) or (cc.nemittx is None)):
            raise ValueError('If sigmax =-1 valid betafx and nemittx MUST be provided!')

        if cc.sigmax == -1:
            gemittx  = cc.nemittx/(beta_rel*gamma_rel)
            cc.sigmax   = np.sqrt(cc.betafx*gemittx+cc.Dx**2*cc.Dp_p**2)
            cc.sigmax   = cc.sigmax.real

        if (cc.sigmay==-1) and ((cc.betafy is None) or (cc.nemitty is None)):
            raise ValueError('If sigmay =-1 valid betafy and nemitty MUST be provided!')

        if cc.sigmay == -1:
            gemitty  = cc.nemitty/(beta_rel*gamma_rel)
            cc.sigmay   = np.sqrt(cc.betafy*gemitty+cc.Dy**2*cc.Dp_p**2)
            cc.sigmay   = cc.sigmay.real

        self.sigmax = cc.sigmax
        self.sigmay = cc.sigmay
        self.x_beam_pos = cc.x_beam_pos
        self.y_beam_pos = cc.y_beam_pos

        self.q_part = cc.q_part
        self.m0_part = cc.m0_part
        self.energy_eV = cc.energy_eV
        self.energy_J = energy_J
        self.beta_rel = beta_rel
        self.Dp_p = cc.Dp_p
        self.nemittx = cc.nemittx
        self.nemitty = cc.nemitty

        self.beam_field_file = cc.beam_field_file

        self.b_spac = cc.b_spac

        self.fact_beam = cc.fact_beam
        self.coast_dens = cc.coast_dens

        self.flag_bunched_beam = cc.flag_bunched_beam

        self.sigmaz = cc.sigmaz
        self.t_offs = cc.t_offs
        self.filling_pattern_file = cc.filling_pattern_file

        self.beam_long_prof_file = cc.beam_long_prof_file

        self.Dh_beam_field = cc.Dh_beam_field
        self.f_telescope_beam = cc.f_telescope_beam
        self.target_grid_beam = cc.target_grid_beam
        self.N_nodes_discard_beam = cc.N_nodes_discard_beam
        self.N_min_Dh_main_beam = cc.N_min_Dh_main_beam

        self.save_beam_field_file_as = cc.save_beam_field_file_as

        self.Nx = cc.Nx
        self.Ny = cc.Ny
        self.nimag = cc.nimag


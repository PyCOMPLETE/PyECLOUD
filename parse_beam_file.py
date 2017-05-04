import numpy as np
from numpy import sqrt

class beam_descr_from_fil:
    def __init__(self, beamfilename, betafx_from_mach_parms_file, Dx_from_mach_parms_file, \
                    betafy_from_mach_parms_file, Dy_from_mach_parms_file):

        sigmax = -1
        sigmay = -1
        sigmaz = -1

        nemittx = None
        nemitty = None


        m0_part =  1.672621777e-27 #kg
        Dp_p = 0.

        Dx = None
        Dy = None
        betafx = None
        betafy = None

        x_beam_pos = 0.
        y_beam_pos = 0.

        Nx = None
        Ny = None
        nimag = None

        Dh_beam_field = None
        f_telescope_beam = None
        target_grid_beam = None
        N_nodes_discard_beam = None
        N_min_Dh_main_beam = None
        save_beam_field_file_as = None

        coast_dens = 0.
        beam_long_prof_file = None

        print 'Parsing beam file: %s'%beamfilename

        f=open(beamfilename)
        exec(f.read())
        f.close()

        # keep beackwards compatibility (optics in the machine params file)
        if Dx == None:
            Dx = Dx_from_mach_parms_file
        if Dy == None:
            Dy = Dy_from_mach_parms_file
        if betafx == None:
            betafx = betafx_from_mach_parms_file
        if betafy == None:
            betafy = betafy_from_mach_parms_file

        qe=1.602176565e-19;
        c=299792458.;

        energy_J = energy_eV * qe
        gamma_rel= energy_J/(m0_part*c*c)
        beta_rel = sqrt(1-1/(gamma_rel*gamma_rel))
        beta_rel = beta_rel.real

        if (sigmax==-1) and ((betafx is None) or (nemittx is None)):
            raise ValueError('If sigmax =-1 valid betafx and nemittx MUST be provided!')

        if sigmax == -1:
            gemittx  = nemittx/(beta_rel*gamma_rel)
            sigmax   = sqrt(betafx*gemittx+Dx*Dx*Dp_p*Dp_p)
            sigmax   = sigmax.real

        if (sigmay==-1) and ((betafy is None) or (nemitty is None)):
            raise ValueError('If sigmay =-1 valid betafy and nemitty MUST be provided!')

        if sigmay == -1:
            gemitty  = nemitty/(beta_rel*gamma_rel)
            sigmay   = sqrt(betafy*gemitty+Dy*Dy*Dp_p*Dp_p)
            sigmay   = sigmay.real

        self.sigmax = sigmax
        self.sigmay = sigmay
        self.x_beam_pos = x_beam_pos
        self.y_beam_pos = y_beam_pos

        self.m0_part = m0_part
        self.energy_eV = energy_eV
        self.energy_J = energy_J
        self.beta_rel = beta_rel
        self.Dp_p = Dp_p
        self.nemittx = nemittx
        self.nemitty = nemitty

        self.beam_field_file = beam_field_file

        self.b_spac = b_spac

        self.fact_beam = fact_beam
        self.coast_dens = coast_dens

        self.flag_bunched_beam = flag_bunched_beam

        self.sigmaz = sigmaz
        self.t_offs = t_offs
        self.filling_pattern_file = filling_pattern_file

        self.beam_long_prof_file = beam_long_prof_file

        self.Dh_beam_field = Dh_beam_field
        self.f_telescope_beam = f_telescope_beam
        self.target_grid_beam = target_grid_beam
        self.N_nodes_discard_beam = N_nodes_discard_beam
        self.N_min_Dh_main_beam = N_min_Dh_main_beam

        self.save_beam_field_file_as = save_beam_field_file_as

        self.Nx = Nx
        self.Ny = Ny
        self.nimag = nimag


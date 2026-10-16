# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #


class Cloud(object):
    def __init__(self, cloudname, config_dict, MP_e, impact_man, dynamics, pyeclsaver,
                 gas_ion_flag, resgasion, t_ion, photoem_flag, phemiss, rho):

        self.name = cloudname
        self.config_dict = config_dict
        self.MP_e = MP_e
        self.impact_man = impact_man
        self.dynamics = dynamics
        self.gas_ion_flag = gas_ion_flag
        self.resgasion = resgasion
        self.t_ion = t_ion
        self.photoem_flag = photoem_flag
        self.phemiss = phemiss
        self.pyeclsaver = pyeclsaver
        self.rho = rho


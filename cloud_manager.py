#----------------------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 7.5.0
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


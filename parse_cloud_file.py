import numpy as np
from scipy.constants import c, e as qe

from . import input_parameters_format_specification as inp_spec
from . import myloadmat_to_obj as mlm


class cloud_descr_from_file:
    def __init__(self, cloudfilename, default_param_obj):

        if cloudfilename is None:
            print('Copying cloud parameters from main input')

            cloudfilename = 'default'
            config_dict = {}
            inp_spec.copy_to_config_dict(config_dict, 'additional_cloud_parameters', default_param_obj)

        else:
            print('Parsing cloud file: %s' % cloudfilename)

            # Parse cloud input file
            cloud_params = inp_spec.import_module_from_file('additional_cloud_parameters', cloudfilename)
        
            # Verify validity of provided module w.r.t. parameters_dict
            inp_spec.assert_module_has_parameters(cloud_params, 'additional_cloud_parameters')
        
            # Create config_dict with cloud parameters (unspecified parameters are set to values in default object)
            config_dict = {}
            inp_spec.update_config_dict(config_dict, cloud_params, 'additional_cloud_parameters', default_obj=default_param_obj)

        config_dict['cloudname'] = cloudfilename 

        cc = mlm.obj_from_dict(config_dict)

        self.cc = cc

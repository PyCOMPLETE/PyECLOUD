from __future__ import division, print_function

import sys
import os
import time
import shutil
import types
from .default_input_parameters import parameters_dict

class PyECLOUD_ConfigException(Exception):
    pass

def assert_module_has_parameters(module, module_name):
    mandatory_parameters = parameters_dict[module_name]['mandatory']
    optional_parameters = parameters_dict[module_name]['optional']
    if isinstance(optional_parameters, dict):
        optional_parameters = set(optional_parameters.keys())

    # Filter out (a) builtin python variables and (b) imported modules
    module_contents = [x for x in dir(module) if (not x.startswith('_'))]
    module_contents = set([x for x in module_contents if (not isinstance(getattr(module, x), types.ModuleType))])

    allowed_parameters = mandatory_parameters | optional_parameters

    missing_parameters = mandatory_parameters - module_contents
    if missing_parameters:
        raise PyECLOUD_ConfigException('Error! These mandatory parameters are not provided by %s: %r' % (module_name, missing_parameters))

    extra_parameters = module_contents - allowed_parameters
    if extra_parameters:
        raise PyECLOUD_ConfigException('Error! These parameters should not be in %s: %r' % (module_name, extra_parameters))


def update_config_dict(config_dict, module, module_name, verbose=True):

    mandatory_parameters = parameters_dict[module_name]['mandatory']
    optional_parameters = parameters_dict[module_name]['optional']


    for parameter in mandatory_parameters:

        value = getattr(module, parameter)
        config_dict[parameter] = value
        if verbose:
            print('%s: %s = %s' % (module_name, parameter, value))

    for parameter, default_value in optional_parameters.items():
        if hasattr(module, parameter):
            if parameter in config_dict and parameter != 't_ion':
                raise PyECLOUD_ConfigException('Parameter %s is specified multiple times!' % parameter)
            value = getattr(module, parameter)
        else:
            if parameter in config_dict:
                if parameter == 't_ion':
                    continue
                else:
                    raise PyECLOUD_ConfigException('Parameter %s is specified multiple times!' % parameter)
            value = default_value
        config_dict[parameter] = value
        if verbose:
            print('%s: %s = %s' % (module_name, parameter, value))


def import_module_from_file(module_name, file_name):
    # Load any file as a python module. This function works for python2 and python3.
    # See https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path
    # "One might think that python imports are getting more and more complicated with each new version."

    # The importlib.util method requires the file name to end with '.py'.
    # Also, loading the module from a temporary directory does not leave any .pyc files.
    dir_name = '/tmp/PyECLOUD_%i' % os.getpid()
    real_module_name = module_name + '_%f' % time.time()
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)

    try:
        new_file_name = dir_name+'/temp_file_%s.py' % (str(time.time()).replace('.','_'))
        shutil.copy(file_name, new_file_name)

        if sys.version_info.major == 2:
            import imp
            module = imp.load_source(real_module_name, new_file_name)

        elif sys.version_info.major == 3:
            import importlib.util
            # In python3, imp.load_source works as well, however it does not return the exact same as the import statement.
            # There are some superfluous contents in the namespace which collides with the inp_spec.assert_module_has_parameters method.

            spec = importlib.util.spec_from_file_location(real_module_name, new_file_name)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
    finally:
        shutil.rmtree(dir_name)

    return module


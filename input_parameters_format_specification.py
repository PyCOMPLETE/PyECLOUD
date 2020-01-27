
import sys
import os
import time
import shutil
import types
from .default_input_parameters import parameters_dict


class PyECLOUD_ConfigException(Exception):
    pass


def assert_module_has_parameters(module, module_name):

    # Format error message
    err_msg = ''
    if module_name == 'combined_simulations_secondaryEmission_machine_parameters':
        nice_module_name = 'simulation_parameters or secondary_emission_parameters or machine_parameters'
    else:
        nice_module_name = module_name

    # Filter out (a) underscorred variables and (b) imported modules
    module_contents = set()
    for attr in dir(module):
        if (not attr.startswith('_')) and (not isinstance(getattr(module, attr), types.ModuleType)):
            module_contents.add(attr)

    # Verify validity of provided module w.r.t. parameters_dict
    mandatory_parameters = parameters_dict[module_name]['mandatory']
    optional_parameters = set(parameters_dict[module_name]['optional'].keys())

    allowed_parameters = set.union(mandatory_parameters, optional_parameters, parameters_dict['superparameters'])

    extra_parameters = set.difference(module_contents, allowed_parameters)
    if extra_parameters:
        err_msg += 'Error! These parameters should not be in %s: %r' % (nice_module_name, extra_parameters)

    missing_parameters = set.difference(mandatory_parameters, module_contents)
    if missing_parameters:
        err_msg += '\nError! These mandatory parameters are not provided by %s: %r' % (nice_module_name, missing_parameters)

    if err_msg:
        raise PyECLOUD_ConfigException(err_msg)


def update_module(module, new_module):
    """
    Similar as the dict.update method, but for modules.
    Some entries are skipped.
    """
    for attr in dir(new_module):

        if attr == 'pi':
            continue

        value = getattr(new_module, attr)
        if attr.startswith('_') or isinstance(value, types.ModuleType):
            continue
        elif hasattr(module, attr):
            raise PyECLOUD_ConfigException('Parameter %s is specified in %s and %s!' % (attr, module.__name__, new_module.__name__))
        else:
            setattr(module, attr, value)


def update_config_dict(config_dict, module, module_name, verbose=False, default_obj=None):

    mandatory_parameters = parameters_dict[module_name]['mandatory']
    optional_parameters = parameters_dict[module_name]['optional']

    for parameter in mandatory_parameters:

        value = getattr(module, parameter)
        config_dict[parameter] = value
        if verbose:
            print('%s: %s = %s' % (module_name, parameter, value))

    for parameter, default_value in list(optional_parameters.items()): # iterates on keys and values of the dictionary

        if hasattr(module, parameter): # the parameter is specified in the input file
            value = getattr(module, parameter)
        elif (default_obj is not None) and hasattr(default_obj, parameter):
            value = getattr(default_obj, parameter)
        else: # the parameter is not specified in the input file (default to be used)
            value = default_value

        # Check for duplicates
        if parameter in config_dict:
            raise PyECLOUD_ConfigException('Parameter %s is specified multiple times!' % parameter)

        # Update config dictionary
        config_dict[parameter] = value

        if verbose:
            print('%s: %s = %r' % (module_name, parameter, value))


def copy_to_config_dict(config_dict, module_name, default_obj, verbose=False):

    mandatory_parameters = parameters_dict[module_name]['mandatory']
    optional_parameters = parameters_dict[module_name]['optional']

    for parameter in mandatory_parameters:
        value = getattr(default_obj, parameter)
        config_dict[parameter] = value

        if verbose:
            print('%s: %s = %s' % (module_name, parameter, value))

    for parameter in optional_parameters:
        value = getattr(default_obj, parameter)
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
    real_module_name = module_name + '_%s' % str(time.time()).replace('.', '_')
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)

    try:
        new_file_name = dir_name + '/temp_file_%s.py' % (module_name + '_' + str(time.time())).replace('.', '_')
        shutil.copy(file_name, new_file_name)

        # As we use pi in many old input files
        with open(new_file_name, 'r') as fid:
            content = fid.read()
        content = "from numpy import pi\n\n" + content
        with open(new_file_name, 'w') as fid:
            fid.write(content)

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

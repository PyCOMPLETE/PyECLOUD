import types
from default_input_parameters import parameters_dict

class PyECLOUD_ConfigException(Exception):
    pass

def assert_module_has_parameters(module, module_name):
    mandatory_parameters = parameters_dict[module_name]['mandatory']
    optional_parameters = parameters_dict[module_name]['optional']
    if isinstance(optional_parameters, dict):
        optional_parameters = set(optional_parameters.keys())

    # Filter out (a) builtin python variables and (b) imported modules
    module_contents = filter(lambda x: (not x.startswith('__')), dir(module))
    module_contents = set(filter(lambda x: (not isinstance(getattr(module, x), types.ModuleType)), module_contents))

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
            print('%s = %s' % (parameter, value))

    for parameter, default_value in optional_parameters.iteritems():
        if hasattr(module, parameter):
            value = getattr(module, parameter)
            if parameter in config_dict and parameter != 't_ion':
                raise PyECLOUD_ConfigException('Parameter %s is specified multiple times!' % parameter)
            config_dict[parameter] = value
        else:
            value = default_value
            if parameter not in config_dict:
                config_dict[parameter] = value
            elif parameter == 't_ion':
                continue
            else:
                raise PyECLOUD_ConfigException('Parameter %s is specified multiple times!' % parameter)

        if verbose:
            print('%s = %s' % (parameter, value))


# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #


def parse_sim_class_string(sim_class_string):
    dict_kwargs = {}
    if '(' in sim_class_string:
        assert(')' in  sim_class_string)

        class_args_string = ')'.join(sim_class_string.split('(', 1)[1].split(')')[:-1])
        module_class_strings = sim_class_string.split('(', 1)[0].split('.')

        module_name = '.'.join(module_class_strings[:-1])
        class_name = module_class_strings[-1]
        for astr in class_args_string.split(','):
            assert('=' in astr)
            nn = astr.split('=')[0].replace(' ', '')
            vv = eval(astr.split('=')[1])
            dict_kwargs[nn] = vv

    else:
        sim_class_strings = sim_class_string.split('.')
        module_name = '.'.join(sim_class_strings[:-1])
        class_name = sim_class_strings[-1]

    return module_name, class_name, dict_kwargs

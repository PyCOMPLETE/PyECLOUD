# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

import importlib
import os

def _load_module(param_file):

    spec = importlib.util.spec_from_file_location('temp.mod',
        param_file)
    mod_params =  importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod_params)
    return mod_params

class SimConfig(object):

    def __init__(self, param_file,
        fname_pkl_at_update='./sim_param.pkl'):

        self.fname_pkl_at_update = fname_pkl_at_update
        self.loaded_files = []

        self.update(param_file)

    def update(self, param_file):

        mod_params = _load_module(param_file)
        abspath_mod = os.path.abspath(param_file)
        absfol_mod = '/'.join(abspath_mod.split('/')[:-1])

        for nn in dir(mod_params):
            if nn.startswith('__') and nn.endswith('__'):
                continue
            setattr(self, nn, getattr(mod_params, nn))

            try:
                vv = getattr(self, nn)
                if '!folder_of_this_file!' in vv:
                    setattr(self, nn,
                        vv.replace('!folder_of_this_file!', absfol_mod))
            except Exception:
                pass

        self.loaded_files.append(os.path.abspath(param_file))

        if self.fname_pkl_at_update is not None:
            self.to_pickle(fname=self.fname_pkl_at_update)

    def to_pickle(self, fname):
        import pickle
        with open(fname, 'wb') as fid:
            pickle.dump(self, fid)

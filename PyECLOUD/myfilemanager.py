import numpy as np


class obj_from_dict:
    def __init__(self, dictto):
        for kk in list(dictto.keys()):
            setattr(self, kk, dictto[kk])


def obj_to_dict(obj):
    dict_out={}
    members = dir(obj)
    for member in members:
        dict_out[member] = getattr(obj, member)
    return dict_out

def myloadmat(filename, squeeze = True):
    import scipy.io as sio
    dict_var=sio.loadmat(filename)
    if squeeze:
        for kk in list(dict_var.keys()):
            try:
                dict_var[kk]=np.squeeze(dict_var[kk])
            except:
                pass
    return dict_var

def myloadmat_to_obj(filename, squeeze = True):
    return  obj_from_dict(myloadmat(filename, squeeze=squeeze))

def dict_of_arrays_and_scalar_from_h5(filename):
    import h5py
    with h5py.File(filename, 'r') as fid:
        f_dict = {}
        for kk in list(fid.keys()):
            f_dict[kk] = np.array(fid[kk]).copy()
            if f_dict[kk].shape == ():
                f_dict[kk] = f_dict[kk].tolist()
    return  f_dict

def object_with_arrays_and_scalar_from_h5(filename):
    return  obj_from_dict(dict_of_arrays_and_scalar_from_h5(filename))


def monitorh5_to_dict(filename, key= 'Bunch'):
    import h5py
    with h5py.File(filename, 'r') as monitor_ev:
        monitor = monitor_ev[key]
        monitor_dict = {}
        for kk in list(monitor.keys()):
            monitor_dict[kk] = np.array(monitor[kk]).copy()

    return monitor_dict


def monitorh5_to_obj(filename, key= 'Bunch'):
    return  obj_from_dict(monitorh5_to_dict(filename, key))

def monitorh5list_to_dict(filename_list, key='Bunch', flag_transpose=False, permissive=False):
    monitor_dict = monitorh5_to_dict(filename_list[0], key=key)
    for i_file in range(1, len(filename_list)):
        print(('Loading '+filename_list[i_file]))
        try:
            monitor_dict_curr = monitorh5_to_dict(filename_list[i_file], key=key)
            if flag_transpose:
                for kk in list(monitor_dict.keys()):
                    monitor_dict[kk] = np.array(list(monitor_dict[kk].T)+list(monitor_dict_curr[kk].T)).T
            else:
                for kk in list(monitor_dict.keys()):
                    monitor_dict[kk] = np.array(list(monitor_dict[kk])+list(monitor_dict_curr[kk]))
        except IOError as err:
            print('Got:')
            print(err)
            if not permissive:
                raise err

    return monitor_dict

def monitorh5list_to_obj(filename_list, key= 'Bunch', flag_transpose=False, permissive=False):
    return  obj_from_dict(monitorh5list_to_dict(filename_list, key, flag_transpose, permissive))


def dict_to_h5(dict_save, filename, compression=None, compression_opts=None):
    import h5py
    with h5py.File(filename, 'w') as fid:
        for kk in list(dict_save.keys()):
            fid.create_dataset(kk, data=dict_save[kk],
		compression=compression, compression_opts=compression_opts)







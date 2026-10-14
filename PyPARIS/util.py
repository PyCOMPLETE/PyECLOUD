# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

from .ring_of_CPUs import RingOfCPUs

class DummyComm(object):
    def __init__(self, N_cores_pretend, pretend_proc_id):
        self.N_cores_pretend = N_cores_pretend
        self.pretend_proc_id = pretend_proc_id

    def Get_size(self):
        return self.N_cores_pretend

    def Get_rank(self):
        return self.pretend_proc_id

    def Barrier(self):
        pass


def get_sim_instance(sim_content, N_cores_pretend,
    id_pretend, init_sim_objects_auto=True):

    myCPUring = RingOfCPUs(
        sim_content,
        comm=DummyComm(N_cores_pretend, id_pretend),
        init_sim_objects_auto=init_sim_objects_auto,
    )
    return myCPUring.sim_content


def get_serial_CPUring(sim_content, init_sim_objects_auto=True):

    myCPUring = RingOfCPUs(
        sim_content, force_serial=True, init_sim_objects_auto=init_sim_objects_auto
    )
    return myCPUring

# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

import multiprocessing as mp
import numpy as np
import os, sys
import importlib
from . import parse_sim_class_string as psc
import traceback
import time

class mpComm(object):
    def __init__(self, pid, N_proc, queue_list,
                    mutex, barriex, turnstile, turnstile2, cnt):
        self._pid = pid
        self._N_proc = N_proc
        self._queue_list = queue_list
        self.mutex = mutex
        self.turnstile = turnstile
        self.turnstile2 = turnstile2
        self.cnt = cnt

    def Get_size(self):
        return self._N_proc

    def Get_rank(self):
        return self._pid

    def Barrier(self):
        self.mutex.acquire()
        self.cnt.value += 1
        if self.cnt.value == self._N_proc:
            self.turnstile2.acquire()
            self.turnstile.release()
        self.mutex.release()
        self.turnstile.acquire()
        self.turnstile.release()
        #criticalpoint
        self.mutex.acquire()
        self.cnt.value -= 1
        if self.cnt.value == 0:
           self.turnstile.acquire()
           self.turnstile2.release()
        self.mutex.release()
        self.turnstile2.acquire()
        self.turnstile2.release()

    def Sendrecv(self, sendbuf, dest, sendtag, recvbuf, source, recvtag):
        self._queue_list[dest].put(sendbuf)
        temp = self._queue_list[self._pid].get()
        recvbuf[:len(temp)]=temp

    def Bcast(self, buf, root):
        self.Barrier()
        if self._pid==root:
            for ii, q in enumerate(self._queue_list):
                if ii!=root:
                    q.put(buf)
        else:
            temp = self._queue_list[self._pid].get()
            buf[:len(temp)]=temp
        self.Barrier()

def todo(sim_module_string, pid, N_proc, queue_list,
                    mutex, barriex, turnstile, turnstile2, cnt, multiturn):
    try:
        comm = mpComm(pid, N_proc, queue_list,
                        mutex, barriex, turnstile, turnstile2, cnt)

        BIN = os.path.expanduser("./")
        sys.path.append(BIN)

        # if len(sim_module_strings)!=2:
        #     raise(ValueError('\n\nsim_class must be given in the form: module.class.\nNested referencing not implemented.\n\n'))

        module_name, class_name, dict_kwargs = psc.parse_sim_class_string(
                sim_module_string)

        SimModule = importlib.import_module(module_name)
        SimClass = getattr(SimModule, class_name)

        simulation_content = SimClass(**dict_kwargs)
        if multiturn:
            from PyPARIS.ring_of_CPUs_multiturn import RingOfCPUs_multiturn
            myCPUring = RingOfCPUs_multiturn(simulation_content, comm=comm)
            myCPUring.run()
        else:
            from PyPARIS.ring_of_CPUs import RingOfCPUs
            myCPUring = RingOfCPUs(simulation_content, comm=comm)
            myCPUring.run()
    except Exception as e:
        print(f"[ERROR] An exception occurred: {e}", file=sys.stderr)
        traceback.print_exc()
        os._exit(1)


if __name__=='__main__':
    import sys

    if len(sys.argv)<4:
        raise ValueError('\n\nSyntax must be:\n\t multiprocexec.py -n N_proc sim_class=module.class\n\n')
    if '-n' not in sys.argv[1] or 'sim_class' not in sys.argv[3]:
        raise ValueError('\n\nSyntax must be:\n\t multiprocexec.py -n N_proc sim_class=module.class\n\n')

    if '--multiturn' in sys.argv:
        multiturn=True
    else:
        multiturn=False

    N_proc = int(sys.argv[2])

    sim_module_string = sys.argv[3].split('=', 1)[1]

    queue_list = [mp.Queue() for _ in range(N_proc)]

    mutex = mp.Semaphore(1)
    barrier = mp.Semaphore(0)
    turnstile = mp.Semaphore(0)
    turnstile2 = mp.Semaphore(1)
    cnt = mp.Value('i', 0)

    proc_list = []
    for pid in range(N_proc):
        proc_list.append(mp.Process(target=todo,
            args=(sim_module_string, pid, N_proc, queue_list,
                    mutex, barrier, turnstile, turnstile2, cnt, multiturn)))
    for p in proc_list:
        p.start()

    accepted_exits = [0, 177] # Sucess exit and resubmit exit
    sim_start_time = time.time()
    # Thread monitor
    while True:
        alive = [p.is_alive() for p in proc_list]
        exitcodes = [p.exitcode for p in proc_list]
        if any(code is not None and code not in accepted_exits for code in exitcodes):
            print("[Main] ERROR: One or more processes failed. Terminating all...", file=sys.stderr)
            for p in proc_list:
                if p.is_alive():
                    p.terminate()
            sys.exit(1)
        if not any(alive):  # all finished
            break
        if time.time()  - sim_start_time > 300: # Stop checks if simulation has ran successfully for more than 5 minutes
            break

        time.sleep(5)  # avoid busy loop, check every 5 seconds
    
    for p in proc_list:
        p.join()
    exitcodes = [p.exitcode for p in proc_list]

    if 177 in exitcodes:
        sys.exit(177)

    if any(code != 0 for code in exitcodes):
        sys.exit(1)



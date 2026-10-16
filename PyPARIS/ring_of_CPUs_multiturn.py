# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

import numpy as np
import time
import sys, os
import socket

from .ring_of_CPUs import SingleCoreComminicator

from . import communication_helpers as ch
from collections import deque

logfilename = 'pyparislog.txt'
def print2logandstdo(message, mode='a+'):
    print(message)
    with open(logfilename, mode) as fid:
        fid.writelines([message+'\n'])


def verbose_mpi_out(message, myid, mpi_verbose, mode='a+'):
    if mpi_verbose:
        t_now = time.mktime(time.localtime())
        time_string = time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(t_now))

        with open('mpi_logfile_cpu%03d.txt'%myid, mode) as fid:
            fid.writelines([time_string + ' - ' + message + '\n'])
        

class RingOfCPUs_multiturn(object):
    def __init__(self, sim_content, N_pieces_per_transfer=1, force_serial = False, comm=None,
                    N_parellel_rings = 1,
                    N_buffer_float_size = 1000000, N_buffer_int_size = 100,
                    verbose = False, mpi_verbose = False, enable_barriers = False,
                    enable_orders_from_master = True):


        self.iteration = -1
        self.sim_content = sim_content
        self.N_turns = sim_content.N_turns
        
        self.N_pieces_per_transfer = N_pieces_per_transfer
        self.N_buffer_float_size = N_buffer_float_size
        self.N_buffer_int_size = N_buffer_int_size
        self.N_parellel_rings = N_parellel_rings
        
        self.verbose = verbose
        self.mpi_verbose = mpi_verbose
        self.enable_barriers = enable_barriers
        
        self.enable_orders_from_master = enable_orders_from_master
        
        if hasattr(sim_content, 'N_pieces_per_transfer'):
            self.N_pieces_per_transfer = sim_content.N_pieces_per_transfer
            
        if hasattr(sim_content, 'N_buffer_float_size'):
            self.N_buffer_float_size = sim_content.N_buffer_float_size
            
        if hasattr(sim_content, 'N_buffer_int_size'):
            self.N_buffer_int_size = sim_content.N_buffer_int_size

        if hasattr(sim_content, 'N_parellel_rings'):
            self.N_parellel_rings = sim_content.N_parellel_rings

        if int(np.mod(self.N_turns, self.N_parellel_rings)) != 0.:
            raise ValueError('Sorry! N_turns needs to be a multiple of N_parellel_rings!')
        
        if hasattr(sim_content, 'verbose'):
            self.verbose = sim_content.verbose
        
        if hasattr(sim_content, 'mpi_verbose'):
            self.mpi_verbose = sim_content.mpi_verbose
        
        if hasattr(sim_content, 'enable_barriers'):
            self.enable_barriers = sim_content.enable_barriers
        
        
        self.sim_content.ring_of_CPUs = self
        
        # choice of the communicator
        if force_serial:
            comm_info = 'Single CPU forced by user.'
            self.comm = SingleCoreComminicator()
            self.comm_type = 'Serial'
        elif comm is not None:
            comm_info = 'Multiprocessing using communicator provided as argument.'
            self.comm = comm
            self.comm_type = 'Multip'
        else:
            comm_info = 'Multiprocessing via MPI.'
            from mpi4py import MPI
            self.comm = MPI.COMM_WORLD
            self.comm_type = 'MPI'
            
        self.verbose_mpi_out = lambda message: verbose_mpi_out(message, self.comm.Get_rank(), 
                                                                   self.mpi_verbose)

        verbose_mpi_out('Debug file (cpu %d)'%(self.comm.Get_rank()), self.comm.Get_rank(), 
                        self.mpi_verbose, mode = 'w')

        self.verbose_mpi_out('Interpreter at %s (cpu %d)'%(sys.executable, self.comm.Get_rank()))
        if self.mpi_verbose and self.comm_type=='MPI':
           import mpi4py
           self.verbose_mpi_out('mpi4py version: %s (cpu %d)'%(mpi4py.__version__, self.comm.Get_rank()))

        self.verbose_mpi_out('Running on %s (cpu %d)'%(socket.gethostname(), self.comm.Get_rank()))         
        
        if self.enable_barriers:
            self.verbose_mpi_out('At barrier 1 (cpu %d)'%self.comm.Get_rank())
            self.comm.Barrier()
            self.verbose_mpi_out('After barrier 1 (cpu %d)'%self.comm.Get_rank())

        #check if there is only one node
        if self.comm.Get_size()==1:
            #in case it is forced by user it will be rebound but there is no harm in that
            self.comm = SingleCoreComminicator()
            
            
        # get info on the grid
        self.N_nodes = self.comm.Get_size()
        self.N_wkrs = self.N_nodes-1
        self.master_id = 0
        self.myid = self.comm.Get_rank()
        self.I_am_a_worker = self.myid!=self.master_id
        self.I_am_the_master = not(self.I_am_a_worker)
        
        # Handle multiturn parallelism
        if np.mod(self.N_nodes, self.N_parellel_rings) !=0:
            raise ValueError('Number of nodes must be a multiple of number of rings!')
            
        self.N_nodes_per_ring = int(float(self.N_nodes)/float(self.N_parellel_rings))
        self.myring = int(np.floor(float(self.myid)/float(self.N_nodes_per_ring)))
        self.myid_in_ring = int(np.mod(float(self.myid),float(self.N_nodes_per_ring)))
        self.I_am_at_start_ring = self.myid_in_ring == 0
        self.I_am_at_end_ring = self.myid_in_ring == (self.N_nodes_per_ring-1)
        
        if self.verbose:
            print2logandstdo("I am %d, master=%s, myring=%d, myid_in_ring=%d (%s%s)"%(
                        self.myid, repr(self.I_am_the_master), self.myring, self.myid_in_ring,
                        {True:'start_ring', False:''}[self.I_am_at_start_ring], 
                        {True:'end_ring', False:''}[self.I_am_at_end_ring]))

        # allocate buffers for communication
        self.buf_float = np.zeros(self.N_buffer_float_size, dtype=np.float64)
        self.buf_int = np.array(self.N_buffer_int_size*[0])

        if self.enable_barriers:
            self.verbose_mpi_out('At barrier 2 (cpu %d)'%self.comm.Get_rank())
            self.comm.Barrier()
            self.verbose_mpi_out('After barrier 2 (cpu %d)'%self.comm.Get_rank())

        if hasattr(self.sim_content, 'pre_init_master'):
            if self.I_am_the_master:
                from_master = self.sim_content.pre_init_master()
            else:
                from_master = None 
            from_master = self._broadcast_from_master(from_master)
            self.sim_content.init_all(from_master)
        else:
            self.sim_content.init_all()
        
        if self.enable_barriers:
            self.verbose_mpi_out('At barrier 3 (cpu %d)'%self.comm.Get_rank())
            self.comm.Barrier()
            self.verbose_mpi_out('After barrier 3 (cpu %d)'%self.comm.Get_rank())

        if self.I_am_the_master:
            print2logandstdo('PyPARIS simulation -- multiturn parallelization')#, mode='w+')
            print2logandstdo(comm_info)
            print2logandstdo('N_cores = %d'%self.N_nodes)
            print2logandstdo('N_pieces_per_transfer = %d'%self.N_pieces_per_transfer)
            print2logandstdo('N_buffer_float_size = %d'%self.N_buffer_float_size)
            print2logandstdo('N_buffer_int_size = %d'%self.N_buffer_int_size)
            print2logandstdo('Multi-ring info:')
            print2logandstdo('N_parellel_rings = %d'%self.N_parellel_rings)
            print2logandstdo('N_nodes_per_ring = %d'%self.N_nodes_per_ring)
            print2logandstdo('Running on %s'%socket.gethostname())
            print2logandstdo('Interpreter at %s'%sys.executable)            
        
        self.left = int(np.mod(self.myid-1, self.N_nodes))
        self.right = int(np.mod(self.myid+1, self.N_nodes))

        if self.enable_barriers:
            self.verbose_mpi_out('At barrier 4 (cpu %d)'%self.comm.Get_rank())
            self.comm.Barrier()
            self.verbose_mpi_out('After barrier 4 (cpu %d)'%self.comm.Get_rank())

        if self.I_am_at_start_ring:
            self.bunches_to_be_treated = deque([])
            self.slices_to_be_treated = []

            self.sim_content.init_start_ring()

        if self.I_am_at_end_ring:
            self.slices_treated = deque([])
        
        if self.I_am_the_master:
            list_bunches = sim_content.init_master()
            self.bunches_to_be_treated.extend(list_bunches)

        if self.enable_barriers:
            self.verbose_mpi_out('At barrier 5 (cpu %d)'%self.comm.Get_rank())
            self.comm.Barrier()
            self.verbose_mpi_out('After barrier 5 (cpu %d)'%self.comm.Get_rank())


    def run(self):
        
        
        self.iteration = 0
        list_received_buffers = [self.sim_content.piece_to_buffer(None)]
        while True:

            orders_from_master = []
            
            ##########################
            # Before slice treatment #
            ##########################
            if self.I_am_at_start_ring:
                # If a bunch is received put in the queue
                buffer_received = list_received_buffers[0]
                bunch_received = self.sim_content.buffer_to_piece(buffer_received)
                if bunch_received is not None:
                    self.bunches_to_be_treated.appendleft(bunch_received)

                # If slices_to_be_treated is empty pop a bunch
                if len(self.slices_to_be_treated)==0 and len(self.bunches_to_be_treated)>0:
                    next_bunch = self.bunches_to_be_treated.pop()
                    
                    # Log some info
                    if self.myring==0 and self.myid_in_ring == 0:
                        t_now = time.mktime(time.localtime())
                        print2logandstdo('%s, iter%03d - cpu %d.%d startin bunch %d/%d turn=%d'%(time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(t_now)), 
                                self.iteration, self.myring, self.myid_in_ring,
                                next_bunch.slice_info['i_bunch'], next_bunch.slice_info['N_bunches_tot_beam'], next_bunch.slice_info['i_turn']))

                
                    self.sim_content.perform_bunch_operations_at_start_ring(next_bunch)

                    if next_bunch.slice_info['i_bunch'] == next_bunch.slice_info['N_bunches_tot_beam']-1:
                        if int(next_bunch.slice_info['i_turn']) >= self.N_turns:
                            orders_from_master.append('stop')
                            if self.I_am_a_worker:
                                raise ValueError('I cannot give orders to anybody :-(')

                    next_bunch.slice_info['i_turn']+=1
                    if next_bunch.slice_info['i_turn'] <= self.N_turns:
                        self.slices_to_be_treated = self.sim_content.slice_bunch_at_start_ring(next_bunch)
                
                # Pop slices
                slice_group = []
                while len(slice_group)<self.N_pieces_per_transfer and len(self.slices_to_be_treated)>0:      
                    slice_group.append(self.slices_to_be_treated.pop())
                if len(slice_group)<1:
                    slice_group.append(None)
            else:
                # Buffer to slice
                slice_group = list(map(self.sim_content.buffer_to_piece, list_received_buffers))
                
            
            
            ####################
            # Treat the slices #
            ####################    
            t_start = time.mktime(time.localtime())              
            for thisslice in slice_group:
                if thisslice is not None:
                    self.sim_content.treat_piece(thisslice)
            t_end = time.mktime(time.localtime()) 
            self._print_some_info_on_comm(thisslice, t_start, t_end)
            
            
            #########################
            # After slice treatment #
            #########################        
            if self.I_am_at_end_ring:    
                
                # Put the slice in slices_treated
                bunch_to_be_sent = None
                for thisslice in slice_group:
                    if thisslice is not None:
                        self.slices_treated.appendleft(thisslice) 
                       
                        if len(self.slices_treated) == self.slices_treated[0].slice_info['N_slices_tot_bunch']:
                            
                            assert(bunch_to_be_sent is None)
                            
                            # Merge slices
                            bunch_to_be_sent = self.sim_content.merge_slices_at_end_ring(self.slices_treated)
        
                            # Perform operations at end ring
                            self.sim_content.perform_bunch_operations_at_end_ring(bunch_to_be_sent)
        
                            # Empty slices_treated
                            self.slices_treated = deque([])
    
                   
                list_of_buffers_to_send = [self.sim_content.piece_to_buffer(bunch_to_be_sent)]
            else:
                # Slice to buffer
                list_of_buffers_to_send = list(map(self.sim_content.piece_to_buffer, slice_group))
            
            
            ########################
            # Send/receive buffer  #
            ########################   
            sendbuf = ch.combine_float_buffers(list_of_buffers_to_send)
            if len(sendbuf) > self.N_buffer_float_size:
                raise ValueError('Float buffer (%d) is too small!\n %d required.'%(self.N_buffer_float_size, len(sendbuf)))

            if self.enable_barriers:
                self.verbose_mpi_out('At barrier L1 (cpu %d)'%self.comm.Get_rank())
                self.comm.Barrier()
                self.verbose_mpi_out('After barrier L1 (cpu %d)'%self.comm.Get_rank())

            self.verbose_mpi_out('At Sendrecv, cpu %d/%d, iter %d'%(self.myid, self.N_nodes, self.iteration))
            self.comm.Sendrecv(sendbuf, dest=self.right, sendtag=self.right, 
                        recvbuf=self.buf_float, source=self.left, recvtag=self.myid)
            self.verbose_mpi_out('After Sendrecv, cpu %d/%d, iter %d'%(self.myid, self.N_nodes, self.iteration))
                
            if self.enable_barriers:
                self.verbose_mpi_out('At barrier L2 (cpu %d)'%self.comm.Get_rank())
                self.comm.Barrier()
                self.verbose_mpi_out('After barrier L2 (cpu %d)'%self.comm.Get_rank())

            list_received_buffers = ch.split_float_buffers(self.buf_float)
            
            
            #######################################################
            # Handle orders (for now only to stop the simulation) #
            #######################################################
            if self.enable_barriers:
                self.verbose_mpi_out('At barrier L3 (cpu %d)'%self.comm.Get_rank())
                self.comm.Barrier()
                self.verbose_mpi_out('After barrier L3 (cpu %d)'%self.comm.Get_rank())


            if self.enable_orders_from_master:
                orders_from_master = self._broadcast_from_master(orders_from_master)
                # check if simulation has to be ended
                if 'stop' in orders_from_master:
                    break

            if self.enable_barriers:
                self.verbose_mpi_out('At barrier L4 (cpu %d)'%self.comm.Get_rank())
                self.comm.Barrier()
                self.verbose_mpi_out('After barrier L4 (cpu %d)'%self.comm.Get_rank())

            self.iteration+=1

            # (TEMPORARY!) To stop
            # if self.iteration==10000:
            #     break
            # (TEMPORARY!)
            
    def _print_some_info_on_comm(self, thisslice, t_start, t_end):
        if self.verbose:
            if thisslice is not None:
                print2logandstdo('Iter start on %s, iter%05d - I am %02d.%02d and I treated slice %d/%d of bunch %d/%d, lasts %ds'%(time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(t_start)), self.iteration, 
                                                self.myring, self.myid_in_ring,
                                                thisslice.slice_info['i_slice'], thisslice.slice_info['N_slices_tot_bunch'], 
                                                thisslice.slice_info['info_parent_bunch']['i_bunch'], 
                                                thisslice.slice_info['info_parent_bunch']['N_bunches_tot_beam'], 
                                                t_end-t_start))
            else:
                print2logandstdo('Iter start on %s, iter%05d - I am %02d.%02d and I treated None, lasts %ds'%(time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(t_start)), self.iteration, 
                                                self.myring, self.myid_in_ring,
                                                t_end-t_start))


    def _broadcast_from_master(self, list_of_strings):
        if self.I_am_the_master:
            # send orders
            buforders = ch.list_of_strings_2_buffer(list_of_strings)
            if len(buforders) > self.N_buffer_int_size:
                raise ValueError('Int buffer is too small!')
            self.buf_int = 0*self.buf_int
            self.buf_int[:len(buforders)] = buforders

            self.verbose_mpi_out('At Bcast, cpu %d/%d, iter %d'%(self.myid, self.N_nodes, self.iteration))
            self.comm.Bcast(self.buf_int, self.master_id)
            self.verbose_mpi_out('After Bcast, cpu %d/%d, iter %d'%(self.myid, self.N_nodes, self.iteration))

        else:    
            # receive orders from the master
            self.verbose_mpi_out('At Bcast, cpu %d/%d, iter %d'%(self.myid, self.N_nodes, self.iteration))
            self.comm.Bcast(self.buf_int, self.master_id)
            self.verbose_mpi_out('After Bcast, cpu %d/%d, iter %d'%(self.myid, self.N_nodes, self.iteration))

        list_of_strings = ch.buffer_2_list_of_strings(self.buf_int)

        return list_of_strings

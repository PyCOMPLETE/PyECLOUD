# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

import h5py

import subprocess
import sys

class SimulationStatus(object):
    def __init__(self,  N_turns_per_run=None, N_turns_target=None, resubmit_command=None, submission_system = None):
        self.N_turns_target = N_turns_target
        self.N_turns_per_run = N_turns_per_run
        self.resubmit_command = resubmit_command
        self.submission_system = submission_system
        
        self.filename = 'simulation_status.sta'
    
    def to_string(self):
        lines = []
        lines.append('present_simulation_part = %d'%self.present_simulation_part)
        lines.append('first_turn_part = %d'%self.first_turn_part)
        lines.append('last_turn_part = %d'%self.last_turn_part)
        lines.append('present_part_done = %s'%repr(self.present_part_done))
        lines.append('present_part_running = %s'%repr(self.present_part_running))
        lines.append('first_run = %s'%repr(self.first_run))

        return '\n'.join(lines)

    def from_string(self, string):
        ddd = {}
        exec(string, ddd)
        self.present_simulation_part = ddd['present_simulation_part']
        self.first_turn_part = ddd['first_turn_part']
        self.last_turn_part = ddd['last_turn_part']
        self.present_part_done = ddd['present_part_done']
        self.present_part_running = ddd['present_part_running']        
        self.first_run = ddd['first_run']

    def dump_to_file(self):
        with open(self.filename, 'w') as fid:
            fid.write(self.to_string())
            
    def load_from_file(self):
        with open(self.filename) as fid:
            self.from_string(fid.read())
        
    def print_from_file(self):
        with open(self.filename) as fid:
            print(fid.read())
    
    def before_simulation(self):
        try:
            self.load_from_file()
            self.present_simulation_part+=1
            self.first_turn_part += self.N_turns_per_run
            self.last_turn_part += self.N_turns_per_run
            self.first_run = False
        except IOError:
            print('Simulation Status not found --> initializing simulation')
            self.present_simulation_part = 0
            self.first_turn_part = 0
            self.last_turn_part = self.N_turns_per_run-1
            self.first_run = True
            self.present_part_done = True
            self.present_part_running = False
            
        if not(self.present_part_done) or self.present_part_running:
            raise ValueError('The previous simulation part seems not finished!!!!')
            
        self.present_part_done = False
        self.present_part_running = True
        
        self.dump_to_file()
        
        print('Starting part:\n\n')
        self.print_from_file()
        print('\n\n')
        
    def after_simulation(self):
        self.load_from_file()
        self.present_part_done = True
        self.present_part_running = False
        self.dump_to_file()
        
        print('Done part:\n\n')
        self.print_from_file()
        print('\n\n')

        if self.last_turn_part+1<self.N_turns_target:
            print('Resubmitting the job.')
            if self.resubmit_command is not None:
                print(f'Resubmit command has been specified: Running {self.resubmit_command}')
                result = subprocess.run(self.resubmit_command, shell=True)
                sys.exit(result.returncode)
            elif self.submission_system is not None:       
                if self.submission_system.lower() == 'slurm':
                    print('Slurm submission selected: Running sbatch 000s_slurm_job')
                    result = subprocess.run(['sbatch', '000s_slurm_job'])
                    sys.exit(result.returncode)
                elif self.submission_system.lower() == 'htcondor':
                    print('HTCondor submission selected: Exiting with resubmit code (177)')
                    sys.exit(177)
                else:
                    raise UserWarning("Unrecognized submission system, cannot resubmit the job")
            else:
                raise UserWarning("No resubmit command or submission system specified, cannot resubmit the job")
                
    def restart_last(self):
        
        self.load_from_file()
        self.N_turns_per_run =  self.last_turn_part - self.first_turn_part + 1

        self.present_simulation_part-=1
        self.first_turn_part -= self.N_turns_per_run
        self.last_turn_part -= self.N_turns_per_run

        self.present_part_done = True
        self.present_part_running = False
        
        self.dump_to_file()
        
        print('Restored status:\n\n')
        self.print_from_file()
        print('\n\n')


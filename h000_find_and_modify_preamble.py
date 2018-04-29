import os

begin =    '''#----------------------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
'''

#-Begin-preamble-------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 7.2.0
#
#
#     Main author:          Giovanni IADAROLA
#                           BE-ABP Group
#                           CERN
#                           CH-1211 GENEVA 23
#                           SWITZERLAND
#                           giovanni.iadarola@cern.ch
#
#     Contributors:         Eleonora Belli
#                           Philipp Dijkstal
#                           Lotta Mether
#                           Annalisa Romano
#                           Giovanni Rumolo
#
#
#     Copyright  CERN,  Geneva  2011  -  Copyright  and  any   other
#     appropriate  legal  protection  of  this  computer program and
#     associated documentation reserved  in  all  countries  of  the
#     world.
#
#     Organizations collaborating with CERN may receive this program
#     and documentation freely and without charge.
#
#     CERN undertakes no obligation  for  the  maintenance  of  this
#     program,  nor responsibility for its correctness,  and accepts
#     no liability whatsoever resulting from its use.
#
#     Program  and documentation are provided solely for the use  of
#     the organization to which they are distributed.
#
#     This program  may  not  be  copied  or  otherwise  distributed
#     without  permission. This message must be retained on this and
#     any other authorized copies.
#
#     The material cannot be sold. CERN should be  given  credit  in
#     all references.
#
#-End-preamble---------------------------------------------------------


for dirpath, _, filenames in os.walk('.'):
    full_paths = [os.path.join(dirpath, x) for x in filenames]

    python_files = [x for x in full_paths if x.endswith('.py')]
    for path in python_files:

        if os.path.abspath(path) == os.path.abspath(__file__):
            continue

        with open(path, 'r') as f:
            content = f.read()
            
        if begin in content:
            print(path)
            content = content.replace(begin, newbegin)
            
            with open(path, 'w') as f:
                f.write(content)
                
        if end in content:
            print 'End to be changed:'
            print(path)
            content = content.replace(end, newend)
            
            with open(path, 'w') as f:
                f.write(content)
                
        if '#--------------' in content:
            print 'Test'
            print(path)
            


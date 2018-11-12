raise ValueError("Not to be used to change the version!\nComment this line to change preamble.")


import os

before_preamble = 0
in_preamble = 1
after_preamble = 2


id_preamble_begin = '#-Begin-preamble----'
id_preamble_end = '#-End-preamble--------'

preamble_new = """#-Begin-preamble-------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
#
#
#     This file is part of the code:
#
#                   PyECLOUD Version 7.6.0
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
#-End-preamble---------------------------------------------------------"""

new_preamble_lines = preamble_new.split('\n')


for dirpath, _, filenames in os.walk('.'):
    full_paths = [os.path.join(dirpath, x) for x in filenames]

    python_files = [x for x in full_paths if x.endswith('.py')]
    for path in python_files:

        if os.path.abspath(path) == os.path.abspath(__file__):
            continue

        with open(path, 'r') as f:
            content = f.read()

        if new_preamble_lines[0] not in content:
            continue

        print(path)

        lines = content.split('\n')

        prev_line = None
        new_lines = []
        status = before_preamble

        for line in lines:
            if status == before_preamble:
                if id_preamble_begin in line:
                    status = in_preamble
                    new_lines.extend(new_preamble_lines)
            elif status == in_preamble:
                if id_preamble_end in prev_line:
                    status = after_preamble
            elif status == after_preamble:
                pass

            if status != in_preamble:
                new_lines.append(line)

            prev_line = line

        new_lines = [x+'\n' for ctr, x in enumerate(new_lines) if ctr+1 != len(new_lines)]

        with open(path, 'w') as f:
            f.writelines(new_lines)
        print('Modified %s' % path)


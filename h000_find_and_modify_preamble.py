import os

begin =    '''#----------------------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
'''

newbegin = '''#-Begin-preamble-------------------------------------------------------
#
#                           CERN
#
#     European Organization for Nuclear Research
'''

end = '''#     The material cannot be sold. CERN should be  given  credit  in
#     all references.
#----------------------------------------------------------------------'''

newend = '''#     The material cannot be sold. CERN should be  given  credit  in
#     all references.
#-End-preamble---------------------------------------------------------'''


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
            


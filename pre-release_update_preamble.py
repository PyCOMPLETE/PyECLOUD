import os
import sys
import re

new_preamble = f'''# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #'''

pattern = re.compile(
    r'#-Begin-preamble-+\n.*?#-End-preamble-+',
    re.DOTALL
)

for dirpath, _, filenames in os.walk('./PyECLOUD'):
    for filename in filenames:
        if not filename.endswith('.py'):
            continue

        path = os.path.join(dirpath, filename)

        if os.path.abspath(path) == os.path.abspath(__file__):
            continue

        with open(path, 'r') as f:
            content = f.read()

        new_content, replacements = pattern.subn(new_preamble, content)

        if replacements:
            print(f'Changing preamble: {path}')

            with open(path, 'w') as f:
                f.write(new_content)
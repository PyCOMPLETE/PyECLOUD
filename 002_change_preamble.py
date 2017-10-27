import sys
import os

before_preamble = 0
in_preamble = 1
after_preamble = 2
in_change = 3

to_add = [
    '#     Author list:          Eleanora Belli',
    '#                           Philipp Dijkstal',
    '#                           Lotta Mether',
    '#                           Annalisa Romano',
]


for dirpath, _, filenames in os.walk('.'):
    full_paths = [os.path.join(dirpath, x) for x in filenames]

    python_files = [x for x in full_paths if x.endswith('.py')]
    for path in python_files:
        with open(path, 'r') as f:
            content = f.read()

        if 'European Organization for Nuclear Research' in content:
            print(path)
        else:
            continue

        lines = content.split('\n')

        prev_line, prev2_line = None, None
        new_lines = []
        status = before_preamble

        for line in lines:
            if status == before_preamble:
                if line.startswith('#-------'):
                    status = in_preamble
            elif status == in_preamble:
                if 'Author and contact' in line:
                    line = line.replace('Author and contact:', 'Main author:       ')
                if 'RUMOLO' in line:
                    status = in_change
                    new_lines.extend(to_add)
            elif status == in_change:
                if 'rumolo@cern.ch' in prev2_line:
                    status = after_preamble
            elif status == after_preamble:
                pass

            if status != in_change:
                new_lines.append(line)

            prev2_line = prev_line
            prev_line = line

        new_lines = [x+'\n' for ctr, x in enumerate(new_lines) if ctr+1 != len(new_lines)]

        #new_file = './test.py'
        with open(path, 'w') as f:
            f.writelines(new_lines)



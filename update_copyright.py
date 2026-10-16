import os
from pathlib import Path

# Bare names match at any depth; paths containing '/' are relative to the
# directory where this script is run. Excluded directories are not traversed.
EXCLUDED_FILES = {
    '_version.py',
    'update_copyright.py',
    '__init__.py',
    '_git.py',
}
EXCLUDED_DIRS = {
    'doc',
    '.git',
    '.venv',
    'venv',
    '__pycache__',
    'testing',
    'docs',
    'other',
    'examples',
    '.github',
    'tests',
}


def is_excluded(path, exclusions):
    relative_path = Path(path).as_posix()
    return any(
        relative_path == Path(entry).as_posix()
        or ('/' not in entry and Path(path).name == entry)
        for entry in exclusions
    )


copyright_statement = (
'''copyright #################################
This file is part of the PyECLOUD Package. 
Copyright (c) CERN, 2011.                  
###########################################''')

config = [
    {'extension': '.py', 'comment_char': '#'},
    {'extension': '.h', 'comment_char': '//'},
    {'extension': '.c', 'comment_char': '//'},
    {'extension': '.pyx', 'comment_char': '#'},
    {'extension': '.f', 'comment_char': '*'},
]

for cc in config:
    extension =  cc['extension']
    comment_char = cc['comment_char']

    cpright_lines = [comment_char + ' ' + line + " " + comment_char + '\n'
                        for line in copyright_statement.splitlines()
                        ] + ['\n']

    for root, dirs, files in os.walk("./"):
        dirs[:] = [
            name for name in dirs
            if not is_excluded(os.path.join(root, name), EXCLUDED_DIRS)
        ]
        for fname in files:
            if fname.endswith(extension):
                file = os.path.join(root, fname)
                if is_excluded(file, EXCLUDED_FILES):
                    continue
                print(file)
                with open(file, 'r') as fid:
                    lines = fid.readlines()
                # Remove copyright statement if present
                if (len(lines) > 0 and
                    lines[0].startswith(comment_char + ' ' + 'copyright ##')):
                    for ill, ll in enumerate(lines):
                        assert ll.startswith(comment_char)
                        if ll.startswith(comment_char + ' ' + '########'):
                            end_cpright = ill
                            break
                    assert lines[end_cpright+1] == '\n'
                    lines = lines[end_cpright+2:]

                lines = cpright_lines + lines

                with open(file, 'w') as fid:
                    fid.writelines(lines)
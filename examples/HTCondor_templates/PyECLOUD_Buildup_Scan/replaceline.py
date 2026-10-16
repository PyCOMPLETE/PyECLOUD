def replaceline_and_save(fname, findln, newline):
    if findln not in newline:
        raise ValueError('Detected inconsistency!!!!')
    
    with open(fname, 'r') as fid:
        lines = fid.readlines()
    
    found = False
    pos = None
    for ii, line in enumerate(lines):
        if line.startswith(findln):
            pos = ii
            found = True
            break
    
    if not found:
        raise ValueError('Not found!!!!')
        
    if '\n' in newline:
        lines[pos] = newline
    else:
        lines[pos] = newline+'\n'
    
    with open(fname, 'w') as fid:
        fid.writelines(lines)

if __name__ == "__main__":
   import argparse
   parser = argparse.ArgumentParser(prog = "replaceline and save")
   parser.add_argument("filename", help="Select path to filename", type = str)
   parser.add_argument("findline", help="Line to replace", type = str)
   parser.add_argument("newline", help="Line to append to file", type = str)
   args = parser.parse_args()
   replaceline_and_save(args.filename,args.findline,args.newline)

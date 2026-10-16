# copyright ################################# #
# This file is part of the PyECLOUD Package.  #
# Copyright (c) CERN, 2011.                   #
# ########################################### #

import os

def htcondor_config(scan_folder, time_requirement_days, job_filename = 'job.job',
        htcondor_files_in='../', runfilename='../run_htcondor', n_cores=1,
	htcondor_subfile='htcondor.sub', listfolderfile = 'list_sim_folders.txt'):

    list_folders= os.listdir(scan_folder)
    list_submit = []
    for folder in list_folders:
        if os.path.isfile(scan_folder + '/' + folder+'/' + job_filename):
            list_submit.append(folder+'\n')
    with open(htcondor_files_in + '/' + listfolderfile, 'w') as fid:
        fid.writelines(list_submit)

    with open(htcondor_files_in + '/' + htcondor_subfile, 'w') as fid:
        fid.write("universe = vanilla\n")
        fid.write("executable = "+ scan_folder+"/$(dirname)/"+job_filename+"\n")
        fid.write('arguments = ""\n')
        fid.write("output = "+ scan_folder+'/$(dirname)/htcondor.out\n')
        fid.write("error = "+scan_folder+"/$(dirname)/htcondor.err\n")
        fid.write("log = "+scan_folder+"/$(dirname)/htcondor.log\n")
        fid.write('transfer_output_files = ""\n')
        fid.write("+MaxRuntime = %d\n"%(time_requirement_days*24*3600))
        fid.write("requestCpus = %d\n"%(n_cores))
        fid.write("queue dirname from %s\n"%listfolderfile)

    with open(htcondor_files_in + '/run_htcondor', 'w') as fid:
        fid.write('condor_submit %s\n'%htcondor_subfile)
        fid.write('condor_q --nobatch\n')
    os.chmod(htcondor_files_in + '/run_htcondor',0o755)

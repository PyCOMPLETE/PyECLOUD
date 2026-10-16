import shutil
import os
import numpy as np
import tarfile

def main(test_folder,default_inputs_folder,folder_prefix="",jobs_fname="jobs.txt"):
    
    os.makedirs(test_folder,exist_ok=True)
    jobs_file = os.path.join(test_folder,jobs_fname)
    if os.path.exists(jobs_file):
        os.remove(jobs_file)
    
    sey_vals = np.arange(1.0,1.61,0.05).round(2)
    intensities = np.arange(0.2,2.31,0.1).round(1)
    b_lengths = np.arange(1.0,1.41,0.1).round(1)
    for b_length in b_lengths:
        for sey_val in sey_vals:
            for intens in intensities:
                generate_config_files(test_folder,default_inputs_folder,sey_val,intens,b_length,folder_prefix,jobs_file)
    

def generate_config_files(test_folder,default_inputs_folder,sey_val,intensity,b_length,folder_prefix,jobs_file):

    config_folder = os.path.join(test_folder,folder_prefix+f"SEY{sey_val:.2f}_intens{intensity:.2f}_blength_{b_length:.2f}")
    os.makedirs(config_folder,exist_ok=True)

    ref_machine_params = os.path.join(default_inputs_folder,"machine_parameters.input")
    config_machine_params = os.path.join(config_folder,"machine_parameters.input")
    shutil.copy(ref_machine_params,config_machine_params)
    

    ref_sim_params = os.path.join(default_inputs_folder,"simulation_parameters.input")
    config_sim_params = os.path.join(config_folder,"simulation_parameters.input")
    shutil.copy(ref_sim_params,config_sim_params)


    ref_sey_params = os.path.join(default_inputs_folder,"secondary_emission_parameters.input")
    config_sey_params = os.path.join(config_folder,"secondary_emission_parameters.input")
    shutil.copy(ref_sey_params,config_sey_params)

    with open(ref_sey_params, "r") as file:
        lines = file.readlines()
    
    with open(config_sey_params, "w") as file:
        for line in lines:
            if "del_max = **" in line:
                file.write(f"del_max = {sey_val:.5f} \n")
            else:
                file.write(line)

    ref_beam_params = os.path.join(default_inputs_folder,"beam.beam")
    config_beam_params = os.path.join(config_folder,"beam.beam")
    shutil.copy(ref_beam_params,config_beam_params)

    with open(config_beam_params, "r") as file:
        lines = file.readlines()
    
    with open(config_beam_params, "w") as file:
        for line in lines:
            if "fact_beam = **" in line:
                file.write(f"fact_beam = {intensity:.5f}e+11 \n")
            elif "sigmaz = **" in line:
                file.write(f"sigmaz = {b_length:.5f}e-09/4.*299792458. \n")
            else:
                file.write(line)
    tar_path = os.path.join(config_folder,"inputs.tgz")

    with tarfile.open(tar_path, "w:gz") as tar:
        for f in os.listdir(config_folder):
            tar.add(os.path.join(config_folder, f), arcname=f)
    
    with open(jobs_file, "a") as file:
        file.write(f"{os.path.relpath(config_folder,'./')}\n")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(prog = "Generate SEY/Intensity/Blength scan")
    parser.add_argument("--prefix",help="Select prefix for sim folder names", default="", type = str)
    parser.add_argument("--jobs_fname",help="Filename containing relpath to all created folders", default="jobs.txt", type = str)
    args = parser.parse_args()
    main("./","./default_input_files",args.prefix,args.jobs_fname)

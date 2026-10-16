#!/bin/bash
set -euo pipefail

# Dict of filling patterns and their full PyECLOUD declaration
declare -A dict
dict["4x72"]="5*(4*(72*[1.]+7*[0.]) + 24*[0.])"
dict["4x36"]="5*(4*(36*[1.]+7*[0.]) + 24*[0.])"

# Configurations to scan
surfaces=("Cu2O" "CuO")
magnets=("dipole" "quadrupole" "drift")
photoems=("pessimistic" "optimistic")

SIM_BASENAME="LHC_6.8TeV_Arcs"
SCRIPT_BASEDIR="$(pwd)"

if [ -f "sims.txt" ]; then
  mv "sims.txt" "last_sim_runs.txt"
fi

# Loop over filling schemes
for key in "${!dict[@]}"; do
	echo "$key -> ${dict[$key]}"
	for surface in "${surfaces[@]}"; do
  		# echo "$surface"
		for magnet in "${magnets[@]}"; do
			# echo "$magnet"
			for photoem in "${photoems[@]}"; do
				# echo "$photoem"
				SIM_BASEFOLDER="${SIM_BASENAME}_${key}_${magnet}_${surface}_${photoem}_photoem"
				echo $SIM_BASEFOLDER
				mkdir -p $SIM_BASEFOLDER/default_input_files
				cp generate_sey_intens_blength_scan.py $SIM_BASEFOLDER

				cp input_templates/beam.beam $SIM_BASEFOLDER/default_input_files/beam.beam
				cp input_templates/simulation_parameters.input $SIM_BASEFOLDER/default_input_files/simulation_parameters.input
				cp input_templates/machine_parameters_${magnet}.input $SIM_BASEFOLDER/default_input_files/machine_parameters.input
				cat input_templates/machine_parameters_${photoem}_photoem.input >> $SIM_BASEFOLDER/default_input_files/machine_parameters.input
				cp input_templates/secondary_emission_parameters_${surface}.input $SIM_BASEFOLDER/default_input_files/secondary_emission_parameters.input
				python replaceline.py $SIM_BASEFOLDER/default_input_files/beam.beam "filling_pattern_file = " "filling_pattern_file = ${dict[$key]}"
				cd $SIM_BASEFOLDER

				python generate_sey_intens_blength_scan.py --prefix $SIM_BASEFOLDER
				cd $SCRIPT_BASEDIR
				cat $SIM_BASEFOLDER/jobs.txt | sed "s/^/${SIM_BASEFOLDER}\//" >> sims.txt
			done
		done
		
	done
done

#!/bin/bash 

set -uo pipefail

ROOT_URL="root://eosproject-e.cern.ch/"
EOS_PATH="/eos/project/e/ecloud-simulations/ekatrali/test_pyparis_resubmit/chroma_00_eldens_1.00e+13/"
SIM_PATH="$EOS_PATH"
CONTAINER_PATH="/cvmfs/unpacked.cern.ch/ghcr.io/ekatralis/ecloud-containers:latest/"
JOB_PWD="$_CONDOR_SCRATCH_DIR"

# Ensure EOS_PATH and SIM_PATH is consistent with trailing slash for later path manipulations
EOS_PATH="${EOS_PATH%/}/"
SIM_PATH="${SIM_PATH%/}/"

xrdcp_opts=(--retry 3)

transfer_inputs() (
	set -euo pipefail

	if xrdfs "$ROOT_URL" stat "${SIM_PATH}simulation_status.sta" >/dev/null 2>&1; then
		xrdcp "${xrdcp_opts[@]}" "${ROOT_URL}${SIM_PATH}simulation_status.sta" ./
		part_num=$(grep -oP 'present_simulation_part\s*=\s*\K\d+' simulation_status.sta)
		part=$(printf "%02d" "$part_num")
		xrdcp "${xrdcp_opts[@]}" "${ROOT_URL}${SIM_PATH}bunch_evolution_${part}.h5" ./
		xrdcp "${xrdcp_opts[@]}" "${ROOT_URL}${SIM_PATH}slice_evolution_${part}.h5" ./
		xrdcp "${xrdcp_opts[@]}" "${ROOT_URL}${SIM_PATH}bunch_status_part${part}.h5" ./
		xrdcp "${xrdcp_opts[@]}" "${ROOT_URL}${SIM_PATH}pyparislog.txt" ./
		if xrdfs "$ROOT_URL" stat "${SIM_PATH}multigrid_config_dip.pkl" >/dev/null 2>&1; then
			xrdcp "${xrdcp_opts[@]}" "${ROOT_URL}${SIM_PATH}multigrid_config_dip.pkl" ./
		fi
		if xrdfs "$ROOT_URL" stat "${SIM_PATH}multigrid_config_dip.txt" >/dev/null 2>&1; then 
			xrdcp "${xrdcp_opts[@]}" "${ROOT_URL}${SIM_PATH}multigrid_config_dip.txt" ./
		fi
		xrdcp "${xrdcp_opts[@]}" "${ROOT_URL}${SIM_PATH}sim_param.pkl" ./
		xrdcp "${xrdcp_opts[@]}" "${ROOT_URL}${SIM_PATH}envinfo.txt" ./
		xrdcp "${xrdcp_opts[@]}" "${ROOT_URL}${SIM_PATH}stdout.txt" ./
		xrdcp "${xrdcp_opts[@]}" "${ROOT_URL}${SIM_PATH}stderr.txt" ./
	else
		echo "No simulation status file found. Initializing simulation..."
	fi

	# LHC_chm_ver.mat: Either here or centrally
	xrdcp "${xrdcp_opts[@]}" "${ROOT_URL}${SIM_PATH}Simulation_parameters.py" ./
	xrdcp "${xrdcp_opts[@]}" -r "${ROOT_URL}${SIM_PATH}pyecloud_config" ./

	apptainer exec --env PYTHONNOUSERSITE=1 --home "$JOB_PWD" "$CONTAINER_PATH" bash -lc 'date; echo $ECLOUD_CONTAINER_VERSION' >> envinfo.txt


)

transfer_outputs() (
	set -euo pipefail
	part_num=$(grep -oP 'present_simulation_part\s*=\s*\K\d+' simulation_status.sta)
	part=$(printf "%02d" "$part_num")
	xrdcp "${xrdcp_opts[@]}" -f "./bunch_evolution_${part}.h5" "${ROOT_URL}${SIM_PATH}"
	xrdcp "${xrdcp_opts[@]}" -f "./slice_evolution_${part}.h5" "${ROOT_URL}${SIM_PATH}"
	xrdcp "${xrdcp_opts[@]}" -f "./bunch_status_part${part}.h5" "${ROOT_URL}${SIM_PATH}"
	xrdcp "${xrdcp_opts[@]}" -f "./pyparislog.txt" "${ROOT_URL}${SIM_PATH}"
	xrdcp "${xrdcp_opts[@]}" -f "./sim_param.pkl" "${ROOT_URL}${SIM_PATH}"
	xrdcp "${xrdcp_opts[@]}" -f "./stdout.txt" "${ROOT_URL}${SIM_PATH}"
	xrdcp "${xrdcp_opts[@]}" -f "./stderr.txt" "${ROOT_URL}${SIM_PATH}"
	if [[ -f "multigrid_config_dip.pkl" ]]; then
		xrdcp "${xrdcp_opts[@]}" -f "./multigrid_config_dip.pkl" "${ROOT_URL}${SIM_PATH}"
	fi
	if [[ -f "multigrid_config_dip.txt" ]]; then
		xrdcp "${xrdcp_opts[@]}" -f "./multigrid_config_dip.txt" "${ROOT_URL}${SIM_PATH}"
	fi
	xrdcp "${xrdcp_opts[@]}" -f "./envinfo.txt" "${ROOT_URL}${SIM_PATH}"
	xrdcp "${xrdcp_opts[@]}" -f "./simulation_status.sta" "${ROOT_URL}${SIM_PATH}"

	
	if (( part_num >= 1 )); then
		prev_part=$(printf "%02d" $((part_num - 1)))
		
		xrdfs ${ROOT_URL} rm \
			"${SIM_PATH}bunch_status_part${prev_part}.h5"
	fi
	
)

transfer_inputs
transfer_status=$?

if (( transfer_status != 0 )); then
    echo "[JOB ERROR]: Input transfer failed"
    exit 1
fi

apptainer exec  --env PYTHONNOUSERSITE=1 \
				--home "$JOB_PWD" \
				"$CONTAINER_PATH" \
				python -m PyPARIS.multiprocexec -n 8 sim_class=PyPARIS_sim_class.Simulation.Simulation \
				>> stdout.txt 2>> stderr.txt
script_exit=$?
echo "Simulation exited with code: $script_exit"

if (( script_exit == 0 || script_exit == 177 )); then
    transfer_outputs
    transfer_status=$?

    if (( transfer_status != 0 )); then
        echo "[JOB ERROR]: Output transfer failed"
        exit 1
    fi
else
    echo "[JOB ERROR]: Simulation failed with exit code ${script_exit}"
fi

exit "$script_exit"

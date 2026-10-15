#!/bin/bash
set -uo pipefail

ROOT_URL="root://eosproject-e.cern.ch/"
SIM_PATH="/eos/project/e/ecloud-simulations/YOURDIRECTORY/stability/run_000/"
# Normalize to one trailing slash before joining URLs and filenames.
while [[ "$ROOT_URL" == */ ]]; do ROOT_URL="${ROOT_URL%/}"; done
ROOT_URL="$ROOT_URL/"
while [[ "$SIM_PATH" == */ ]]; do SIM_PATH="${SIM_PATH%/}"; done
SIM_PATH="$SIM_PATH/"
previous_bunch_state=""
JOB_PWD="${_CONDOR_SCRATCH_DIR:?}"
cpus=${1:?Missing CPU count}
xrdcp_opts=(--retry 3)
cd "$JOB_PWD" || exit 1

run_python() {
    python "$@"
}

checkpoint_files() {
    local part_num part
    part_num=$(sed -n 's/^present_simulation_part = \([0-9][0-9]*\)$/\1/p' simulation_status.sta) || return 1
    [[ "$part_num" =~ ^[0-9]+$ ]] || return 1
    part=$(printf '%02d' "$((10#$part_num))") || return 1
    files=("bunch_evolution_${part}.h5" "slice_evolution_${part}.h5"
           "bunch_status_part${part}.h5" pyparislog.txt sim_param.pkl
           envinfo.txt stdout.txt stderr.txt)
}

transfer_inputs() {
    local listing file
    # A failed listing is an error, not evidence of a fresh simulation.
    listing=$(xrdfs "$ROOT_URL" ls "${SIM_PATH%/}") || return 1
    if grep -Fxq "${SIM_PATH}simulation_status.sta" <<< "$listing"; then
        xrdcp "${xrdcp_opts[@]}" "${ROOT_URL}${SIM_PATH}simulation_status.sta" ./ || return 1
        checkpoint_files || return 1
        for file in "${files[@]}"; do
            xrdcp "${xrdcp_opts[@]}" "${ROOT_URL}${SIM_PATH}${file}" ./ || return 1
        done
        previous_bunch_state="${files[2]}"
        for file in multigrid_config_dip.pkl multigrid_config_dip.txt; do
            if grep -Fxq "${SIM_PATH}${file}" <<< "$listing"; then
                xrdcp "${xrdcp_opts[@]}" "${ROOT_URL}${SIM_PATH}${file}" ./ || return 1
            fi
        done
    else
        echo "No simulation status in the accessible run directory; starting fresh."
    fi
    xrdcp "${xrdcp_opts[@]}" "${ROOT_URL}${SIM_PATH}Simulation_parameters.py" ./ || return 1
    xrdcp "${xrdcp_opts[@]}" -r "${ROOT_URL}${SIM_PATH}pyecloud_config" ./ || return 1
}

transfer_outputs() {
    local file
    grep -Fxq 'present_part_done = True' simulation_status.sta || return 1
    grep -Fxq 'present_part_running = False' simulation_status.sta || return 1
    checkpoint_files || return 1
    # Check all required files before starting publication.
    for file in "${files[@]}"; do
        [[ -f "$file" ]] || { echo "Missing checkpoint file: $file" >&2; return 1; }
    done
    for file in "${files[@]}"; do
        xrdcp "${xrdcp_opts[@]}" -f "$file" "${ROOT_URL}${SIM_PATH}${file}" || return 1
    done
    for file in multigrid_config_dip.pkl multigrid_config_dip.txt; do
        if [[ -f "$file" ]]; then
            xrdcp "${xrdcp_opts[@]}" -f "$file" "${ROOT_URL}${SIM_PATH}${file}" || return 1
        fi
    done
    # Publish the status only after every associated file has arrived.
    xrdcp "${xrdcp_opts[@]}" -f simulation_status.sta \
        "${ROOT_URL}${SIM_PATH}simulation_status.sta" || return 1

    # Delete only the state we resumed from, after the new checkpoint is saved.
    if [[ -n "$previous_bunch_state" && "$previous_bunch_state" != "${files[2]}" ]]; then
        if ! xrdfs "$ROOT_URL" rm "${SIM_PATH}${previous_bunch_state}"; then
            echo "Warning: checkpoint saved, but could not remove $previous_bunch_state" >&2
        fi
    fi
    return 0
}

transfer_inputs || exit 1
# Remove if using container -------------------------------------------------
# Batch shells do not run 'conda init'; load Conda's shell functions first.
eval "$(conda shell.bash hook)" || exit 1
conda activate /eos/user/{yourinitial}/{yourusername}/pyecloud-env || exit 1
# ---------------------------------------------------------------------------
# Append one record per simulation part: start time and software version.
date >> envinfo.txt
run_python -c 'import os, PyECLOUD; key = "ECLOUD_CONTAINER_VERSION"; print(f"{key}={os.environ[key]}" if key in os.environ else f"PyECLOUD={PyECLOUD.__version__}")' >> envinfo.txt

run_python -m PyPARIS.multiprocexec -n "$cpus" \
    sim_class=PyPARIS_sim_class.Simulation.Simulation >> stdout.txt 2>> stderr.txt
script_exit=$?
echo "Simulation exited with code: $script_exit"
if (( script_exit == 0 || script_exit == 177 )); then
    transfer_outputs || exit 1
else
    echo "Simulation failed; completed checkpoint not published." >&2
    tail -n 100 stderr.txt >&2
fi
exit "$script_exit"

#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

mkdir -p outputs
failed_file="outputs/failed_simulations.txt"
: > "$failed_file"
failed=0

cleanup_generated() {
    find "$script_dir" -maxdepth 1 -type f \
        \( -name 'simulation_status.sta' -o -name '*.pkl' -o -name '*.h5' -o -name '*.txt' \) \
        -print -delete
}

run_and_check() {
    local simulation_script="$1"
    local output_dir="outputs/${simulation_script%.sh}"
    echo "Running $simulation_script"

    if bash "$simulation_script"; then
        if ! python 001_check_output.py --output-dir "$output_dir"; then
            echo "Output check failed after $simulation_script" | tee -a "$failed_file"
            failed=1
        fi
    else
        echo "Simulation failed: $simulation_script (exit $?)" | tee -a "$failed_file"
        failed=1
    fi

    echo "Cleaning generated simulation files"
    cleanup_generated
}

# Remove stale simulation files before starting, then run each mode in turn.
cleanup_generated
run_and_check 000a_serial_run_sim_with_eclouds.sh
run_and_check 000b_multiproc_run_sim_with_eclouds.sh
run_and_check 000c_mpi_run_sim_with_eclouds.sh

if [ "$failed" -ne 0 ]; then
    echo "One or more runs failed; see $failed_file"
    exit 1
fi

echo "All simulations and output checks completed successfully."

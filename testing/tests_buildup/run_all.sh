#!/usr/bin/env bash
set -u

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

failed_file="$script_dir/failed_tests.txt"
: > "$failed_file"
failed=0

run_test() {
    local label="$1"
    shift
    echo "Running $label"
    if ! "$@"; then
        printf '%s\n' "$label failed" | tee -a "$failed_file"
        failed=1
    fi
}

run_test "buildup simulations and reference comparisons" python 002_run_all_tests.py
# 002_run_all_tests.py records individual case failures but currently exits successfully.
if [ -s "$failed_file" ]; then
    failed=1
fi

bfield_dir="$script_dir/Bfield_from_file"
run_test "B-field map generation" bash -c 'cd "$1" && python 000_make_bfile.py' _ "$bfield_dir"
run_test "B-field map loading" bash -c 'cd "$1" && python 001_load_bfile.py' _ "$bfield_dir"
run_test "Boris pusher with B-field map" bash -c 'cd "$1" && python 002_test_pusher.py' _ "$bfield_dir"

run_test "uniform vs non-uniform comparison" python 004_compare_unif_vs_non_unif.py

if [ "$failed" -ne 0 ]; then
    echo "One or more tests failed; see $failed_file"
    exit 1
fi

echo "All buildup tests completed successfully."

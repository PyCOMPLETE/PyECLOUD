#!/usr/bin/env bash
set -u

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

mkdir -p outputs
failed_tests="outputs/failed_tests.txt"
: > "$failed_tests"

for test in \
    000_ecloud_instab_sim \
    001_compare_kick_against_headtail \
    002_compare_kick_and_betatron_motion \
    003_particle_tune_shift_with_frozen_cloud \
    004_particle_tune_shift_with_real_cloud \
    005a_particle_tune_shift_comparison_gen_HT_file \
    005b_particle_tune_shift_against_HT \
    005c_particle_tune_shift_against_HT_slice_mode \
    006_test_slice_by_slice_mode \
    007_time_slice_by_slice_mode \
    008_test_multigrid_pinch \
    009_particle_tune_shift_against_HT_multigrid \
    010_test_pyecl_output \
    011_test_multigrid_pinch_record_grid
do
    mkdir -p "outputs/$test"
    echo "Running $test.py"
    if ! python "$test.py"; then
        printf '%s\n' "$test.py" >> "$failed_tests"
    fi
done

if [ -s "$failed_tests" ]; then
    echo "Some tests failed. See $failed_tests"
    exit 1
fi

echo "All tests completed successfully."

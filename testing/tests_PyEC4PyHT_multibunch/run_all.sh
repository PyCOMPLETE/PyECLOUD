#!/usr/bin/env bash
set -u

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
mkdir -p outputs
failed_tests="outputs/failed_tests.txt"
: > "$failed_tests"

for test in \
    000_test_cloudsim \
    001_comparison_against_reference
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

#!/usr/bin/env bash
set -u

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

failed_tests="failed_tests.txt"
: > "$failed_tests"
rm -rf outputs
rm -f output.png

if ! mpiexec -n 4 python simulation_check.py; then
    printf '%s\n' "mpiexec -n 4 python simulation_check.py failed" >> "$failed_tests"
fi

mkdir -p outputs
if [ -f output.png ]; then
    if ! mv output.png outputs/output.png; then
        printf '%s\n' "Could not move output.png into outputs/" >> "$failed_tests"
    fi
else
    printf '%s\n' "simulation_check.py did not create output.png" >> "$failed_tests"
fi

if [ -s "$failed_tests" ]; then
    echo "Parallel buildup test failed; see $failed_tests"
    exit 1
fi

echo "Parallel buildup test completed successfully."

#!/usr/bin/env bash
set -u

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
mkdir -p outputs

if [ ! -f CLIC_DR_n260_optics.h5 ]; then
    python -c 'import pickle; import PyECLOUD.myfilemanager as mfm; optics = pickle.load(open("CLIC_DR_n260_optics.pkl", "rb")); optics.pop("circumference", None); mfm.dict_to_h5(optics, "CLIC_DR_n260_optics.h5")'
fi

failed_tests="outputs/failed_tests.txt"
: > "$failed_tests"

for test in \
    000a_ion_simulation \
    000b_plot_bunch_evolution \
    001a_ion_simulation_multiturn \
    001b_plot_turn_evolution
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

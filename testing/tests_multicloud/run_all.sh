#!/usr/bin/env bash

rm -f failed_tests.txt

for script in *.py; do
    python "$script"
    if [ $? -ne 0 ]; then
        echo "python $script" >> failed_tests.txt
    fi
done

mkdir -p cross_sections/comparison_plots
if ! mv Cross*.png cross_sections/comparison_plots/; then
    echo "Could not move Cross*.png into cross_sections/comparison_plots/" >> failed_tests.txt
fi

if [ -s failed_tests.txt ]; then
    echo "One or more tests failed; see failed_tests.txt"
    exit 1
fi

echo "All tests completed successfully."

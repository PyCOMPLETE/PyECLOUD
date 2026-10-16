#!/usr/bin/env bash

#rm -f failed_tests.txt
rm -r  000_test_photoemission

for test in \
    001_test_multipole \
    002_test_Boris \
    003_test_Boris_quad_compar \
    004_test_Boris_quad \
    005_test_parab_low_ener \
    006_test_ECLOUD_energy_distributions \
    007_test_ECLOUD_sey_curves \
    008_test_ECLOUD_emission_angles \
    009_test_Furman_Pivi_energy_distributions \
    010_test_Furman_Pivi_sey_curves \
    010b_test_Furman_Pivi_sey_curves \
    011_test_Furman_Pivi_model \
    012_test_Furman_Pivi_emission_angles
do
    rm -r "$test"
    # python "$test.py" -o "$test/output" --noshow
    if [ $? -ne 0 ]; then
            echo "python $test.py -o "$test/output" --noshow" >> failed_tests.txt
    fi
done

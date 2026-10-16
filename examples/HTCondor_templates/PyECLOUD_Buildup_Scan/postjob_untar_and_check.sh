#!/usr/bin/env bash
set -u

BASE_DIR="/eos/project/e/ecloud-simulations/ekatrali/LHC_6.8TeV_Arcs_Filling_Pattern_Scan"
SIM_LIST="sims.txt"
FAILED_LIST="failed_sims.txt"

> "$FAILED_LIST"

TOTAL=$(grep -cve '^\s*$' "$SIM_LIST")
COUNT=0

python -c "import os; import scipy.io"
if [[ $? -ne 0 ]]; then
    echo "Error: Python with scipy is required to run this script."
    exit 1
fi

while IFS= read -r relpath; do
    [[ -z "$relpath" ]] && continue

    simdir="$BASE_DIR/$relpath"
    archive="$simdir/output.tgz"
    output_mat="$simdir/output.mat"
    progress_txt="$simdir/progress.txt"

    failed=0
    
    ((COUNT++))
    printf "\r[%d/%d] Processing: %s" "$COUNT" "$TOTAL" "$relpath"

    if [[ ! -d "$simdir" ]]; then
        echo "$simdir" >> "$FAILED_LIST"
        continue
    fi

    if [[ ! -f "$archive" ]]; then
        echo "$simdir" >> "$FAILED_LIST"
        continue
    fi

    tar -xzf "$archive" -C "$simdir"
    if [[ $? -ne 0 ]]; then
        echo "$simdir" >> "$FAILED_LIST"
        continue
    fi

    if [[ ! -f "$output_mat" ]]; then
        failed=1
    fi

    if [[ ! -f "$progress_txt" ]]; then
        failed=1
    else
        progress=$(cat "$progress_txt" | tr -d '[:space:]')

        if ! awk -v p="$progress" 'BEGIN { exit !(p ~ /^[0-9]+(\.[0-9]+)?$/ && p > 0.9) }'; then
            failed=1
        fi
    fi

    python -c "import os; import scipy.io; mat = scipy.io.loadmat('$output_mat'); print(mat.keys())" > /dev/null 2>&1
    if [[ $? -ne 0 ]]; then
        failed=1
    fi

    if [[ "$failed" -eq 1 ]]; then
        echo "$simdir" >> "$FAILED_LIST"
    else
        rm -f "$archive"
    fi

done < "$SIM_LIST"

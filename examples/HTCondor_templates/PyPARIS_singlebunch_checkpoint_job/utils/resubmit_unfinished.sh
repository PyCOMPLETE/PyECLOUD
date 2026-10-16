#!/usr/bin/env bash

set -uo pipefail

BASE_DIR=$(pwd)
CONDOR_CURRENT_STATUS_FILE=".$(date +%Y%m%d)_condor_out"

touch "$CONDOR_CURRENT_STATUS_FILE"
condor_q > "$CONDOR_CURRENT_STATUS_FILE"

for chroma in $(seq -w 5 5 10); do
    for job in $(ls chroma_${chroma}/sims/); do
        # echo "$job"
        job_dir=${BASE_DIR}/chroma_${chroma}/sims/${job}
        if [[ -d "$job_dir" ]]; then
            tail -n 20 "$job_dir/stdout.txt" | grep "Stop simulation due to beam losses." > /dev/null
            if [[ $? -eq 1 ]]; then
		grep -o "${job}" "${CONDOR_CURRENT_STATUS_FILE}" > /dev/null
		if [[ $? -eq 1 ]]; then
                	echo "cd /afs/cern.ch/work/e/ekatrali/private/submit_mirror${job_dir} && mkdir old_submissions && mv ${job}.dag.* old_submissions/ && mv log.* old_submissions/ && condor_submit_dag -batch-name symchroma_${job} ${job}.dag"
            	fi
	    fi
        fi
    done
done

rm "$CONDOR_CURRENT_STATUS_FILE"

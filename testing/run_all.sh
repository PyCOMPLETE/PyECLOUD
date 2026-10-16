#!/usr/bin/env bash
set -u

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
failed_file="$script_dir/run_all_failures.txt"
: > "$failed_file"
failed=0

mapfile -t runners < <(find "$script_dir" -mindepth 2 -maxdepth 2 -type f -name run_all.sh -print | sort)

if [ "${#runners[@]}" -eq 0 ]; then
    echo "No subdirectory run_all.sh scripts found under $script_dir"
    exit 1
fi

for runner in "${runners[@]}"; do
    relative_runner="${runner#"$script_dir"/}"
    runner_dir="$(dirname -- "$runner")"
    echo "Running $relative_runner"
    if (cd -- "$runner_dir" && bash ./run_all.sh); then
        echo "Completed $relative_runner"
    else
        status=$?
        printf '%s (exit %s)\n' "$relative_runner" "$status" | tee -a "$failed_file"
        failed=1
    fi
done

if [ "$failed" -ne 0 ]; then
    echo "One or more subdirectory runners failed; see $failed_file"
    exit 1
fi

echo "All subdirectory runners completed successfully."

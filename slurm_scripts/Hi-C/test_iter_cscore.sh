#!/bin/sh
#SBATCH --job-name=cscore
#SBATCH --mem=50gb
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=cscore_%A-%a.log

# ========== VARIABLES ==========
describer=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)
array_id=${SLURM_ARRAY_TASK_ID}
source ./config.sh

N_RUNS=10

for dir in "${path_cscore}" ; do
  if [ ! -d "${dir}" ]; then
    mkdir -p "${dir}"
  fi
done

# ========== MODULES ==========
module load cscoretool/1.1

echo " ................................................................ START cscore ${describer} ................................................................ "

best_L=-999999999999
best_run=-1

for run in $(seq 1 ${N_RUNS}); do

    echo "  --- Run ${run}/${N_RUNS} ---"

    run_dir="${path_cscore}/${describer}_run${run}"
    mkdir -p "${run_dir}"
    log_file="${run_dir}/cscore_run${run}.log"

    CscoreTool1.1 ${genome100kb} \
        homertxt_${array_id}.txt \
        ${run_dir}/${describer} \
        4 1000000 2>&1 | tee "${log_file}"

    # Extract last L value from log
    last_L=$(grep '^L=' "${log_file}" | tail -1 | sed 's/L=\([0-9.]*\).*/\1/')

    echo "  Run ${run} final L = ${last_L}"

    # Compare with best using awk for float comparison
    is_better=$(awk -v new="${last_L}" -v best="${best_L}" 'BEGIN {print (new > best) ? 1 : 0}')

    if [ "${is_better}" -eq 1 ]; then
        best_L=${last_L}
        best_run=${run}
    fi

done

echo "  Best run: ${best_run} with L = ${best_L}"

# Rename best run folder to final output directory
mv "${path_cscore}/${describer}_run${best_run}" "${path_cscore}/${describer}"
echo "  Best run (run ${best_run}) renamed to ${final_dir}"

# Remove other runs
for run in $(seq 1 ${N_RUNS}); do
    if [ "${run}" -ne "${best_run}" ]; then
        rm -rf "${path_cscore}/${describer}_run${run}"
    fi
done
echo "  Removed suboptimal runs"

echo " ................................................................ END cscore ${describer} ................................................................ "

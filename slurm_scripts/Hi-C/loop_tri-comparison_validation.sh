#!/bin/bash
#SBATCH --job-name=coolpup_comp
#SBATCH --mem=40gb
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --array=1-7
#SBATCH --output=coolpup_comp_%A_%a.log

module load coolpuppy/1.1.0-foss-2021b
module load PyTables/3.6.1-foss-2021b

#############################################
# PATHS
#############################################
BASE="/ijc/LABS/STIK/RAW/Hi-C/HiC_marta_Oct25"
INTERSECT="${BASE}/loops/intersected_loops"
COOL="${BASE}/coolMatrix"
OUT="${INTERSECT}/pup_results/comparisons"
mkdir -p "${OUT}"

#############################################
# PARAMETERS
#############################################
RES="10kb"
COOL_SUFFIX="_${RES}_KR.cool"

#############################################
# SAMPLE DEFINITIONS
#############################################
declare -A COOL_FILES=(
    ["EP1"]="HiC-REH-EP1"
    ["Aro"]="HiC-REH-EP1-Aro"
    ["REH"]="REH"
)

# Always plot in all three conditions
ALL_SAMPLES=("EP1" "Aro" "REH")

#############################################
# GET BEDPE FILE FOR THIS TASK
#############################################
bedpe=$(ls "${INTERSECT}"/*_"${RES}"_*.bedpe | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [[ ! -f "${bedpe}" ]]; then
    echo "[ERROR] BED file not found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

bedpe_name=$(basename "${bedpe}")
comp="${bedpe_name%.bedpe}"
comp="${comp%%_05_*}"

echo "[INFO] =========================================="
echo "[INFO] Task ${SLURM_ARRAY_TASK_ID}: ${comp}"
echo "[INFO] BEDPE file: ${bedpe_name}"
echo "[INFO] Will plot in all samples: ${ALL_SAMPLES[*]}"
echo "[INFO] =========================================="

#############################################
# RUN COOLPUP FOR ALL THREE SAMPLES
#############################################
for sample in "${ALL_SAMPLES[@]}"
do
    cool_name="${COOL_FILES[$sample]}"
    cool_file="${COOL}/${cool_name}${COOL_SUFFIX}"
    out_clpy="${OUT}/${comp}_in_${sample}.clpy"
    out_pdf="${OUT}/pup_${comp}_in_${sample}.pdf"

    # Skip if output already exists
    if [[ -f "${out_clpy}" && -f "${out_pdf}" ]]; then
        echo "[INFO] Output already exists for ${comp} in ${sample}, skipping"
        continue
    fi

    if [[ ! -f "${cool_file}" ]]; then
        echo "[ERROR] Missing matrix for ${sample} (${cool_name})"
        echo "[ERROR] Expected file: ${cool_file}"
        exit 1
    fi

    echo "[INFO] Running coolpup: ${comp} in ${sample}"

    # coolpup.py \
    #     "${cool_file}" \
    #     "${bedpe}" \
    #     --format bedpe \
    #     -o "${out_clpy}"

    # if [[ $? -ne 0 ]]; then
    #     echo "[ERROR] coolpup.py failed for ${comp} in ${sample}"
    #     exit 1
    # fi

    plotpup.py \
        --input_pups "${out_clpy}" \
        --font_scale 0.5 \
        --output "${out_pdf}" \
        --vmin 0.5 \
        --vmax 4.6

    if [[ $? -ne 0 ]]; then
        echo "[ERROR] plotpup.py failed for ${comp} in ${sample}"
        exit 1
    fi

    echo "[INFO] ✓ Completed: ${comp} in ${sample}"
done

echo "[INFO] =========================================="
echo "[INFO] DONE: ${comp}"
echo "[INFO] =========================================="

#!/bin/bash
#SBATCH --job-name=coolpup_cond
#SBATCH --mem=40gb
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=coolpup_cond_%j.log
module load coolpuppy/1.1.0-foss-2021b
module load PyTables/3.6.1-foss-2021b
#############################################
# PATHS
#############################################
BASE="/ijc/LABS/STIK/RAW/Hi-C/HiC_marta_Oct25"
LOOPS="${BASE}/loops"
COOL="${BASE}/coolMatrix"
OUT="${LOOPS}/pup_results/per_condition"
CROSS="${LOOPS}/pup_results/cross_condition"   # <-- new output dir for cross comparisons
TMP="${LOOPS}/pup_results/tmp"
mkdir -p "${OUT}" "${CROSS}" "${TMP}"
#############################################
# PARAMETERS
#############################################
RES="10kb"
COOL_SUFFIX="_${RES}_KR.cool"
#############################################
# VALIDATE LOOPS PER CONDITION
#############################################
echo "========== VALIDATING LOOPS PER CONDITION =========="
for tsv in ${LOOPS}/hic_*_loops_05_${RES}.tsv
do
    sample=$(basename "${tsv}" | sed 's/^hic_//' | sed "s/_loops_05_${RES}.tsv//")
    cool_file="${COOL}/${sample}${COOL_SUFFIX}"
    if [[ ! -f "${cool_file}" ]]; then
        echo "[WARN] Missing matrix for ${sample}, skipping"
        continue
    fi
    bedpe="${TMP}/${sample}_loops_05_${RES}.bedpe"
    # Convert TSV → BEDPE (drop header)
    sed 1d "${tsv}" > "${bedpe}"
    echo "[INFO] coolpup: ${sample} vs its own loops"
#     coolpup.py \
#         "${cool_file}" \
#         "${bedpe}" \
#         -o "${OUT}/${sample}_loops_05_${RES}.clpy"
    plotpup.py \
        --input_pups "${OUT}/${sample}_loops_05_${RES}.clpy" \
        --font_scale 0.5 \
        --vmin 0.5 \
        --vmax 4.6 \
        --output "${OUT}/pup_${sample}_loops_05_${RES}.pdf"
done
echo "========== DONE (CONDITIONS) =========="

#############################################
# CROSS-CONDITION: loops of A applied to matrix of B (A != B)
#############################################
echo "========== CROSS-CONDITION PUP =========="
for tsv in ${LOOPS}/hic_*_loops_05_${RES}.tsv
do
    loop_sample=$(basename "${tsv}" | sed 's/^hic_//' | sed "s/_loops_05_${RES}.tsv//")
    bedpe="${TMP}/${loop_sample}_loops_05_${RES}.bedpe"

    # Ensure bedpe exists (may not if matrix was missing above, so regenerate safely)
    if [[ ! -f "${bedpe}" ]]; then
        sed 1d "${tsv}" > "${bedpe}"
    fi

    for cool_file in ${COOL}/*${COOL_SUFFIX}
    do
        matrix_sample=$(basename "${cool_file}" | sed "s/${COOL_SUFFIX}//")

        # Skip same-sample pair (already done above)
        if [[ "${loop_sample}" == "${matrix_sample}" ]]; then
            continue
        fi

        tag="loops_${loop_sample}__matrix_${matrix_sample}"
        clpy_out="${CROSS}/${tag}_${RES}.clpy"
        pdf_out="${CROSS}/pup_${tag}_${RES}.pdf"

        echo "[INFO] coolpup: loops from ${loop_sample} → matrix ${matrix_sample}"
        # coolpup.py \
        #     "${cool_file}" \
        #     "${bedpe}" \
        #     -o "${clpy_out}"

        plotpup.py \
            --input_pups "${clpy_out}" \
            --font_scale 0.5 \
            --vmin 0.5 \
            --vmax 4.6 \
            --output "${pdf_out}"
    done
done
echo "========== DONE (CROSS-CONDITION) =========="

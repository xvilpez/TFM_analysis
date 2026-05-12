#!/bin/bash
#SBATCH --job-name=loop_chip_coverage
#SBATCH --output=loop_chip_coverage_%j.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00

# ============================================================================
# Loop anchor x ChIP-seq coverage — bedtools
# E2A-PBX1 ChIP on REH/EP1/ARO loop anchors
# ============================================================================

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
ANCHORS_DIR="/ijc/LABS/STIK/RAW/Hi-C/HiC_marta_Oct25/loops/intersected_loops"
BAM="/ijc/LABS/STIK/DATA/CHIP-seq/public_data_september_2023/bowtie2_results/final_data/697HF_E2APBX1-HF_clean.bam"
OUT_DIR="${ANCHORS_DIR}/chip_coverage"
TMP_DIR="${OUT_DIR}/tmp"

mkdir -p "${OUT_DIR}" "${TMP_DIR}"

# ----------------------------------------------------------------------------
# Modules
# ----------------------------------------------------------------------------
module load BEDTools/2.31.1-GCC-11.2.0
module load SAMtools/1.19.2-foss-2021b

# ----------------------------------------------------------------------------
# 1. Verify BAM and index
# ----------------------------------------------------------------------------
echo "[$(date)] Verifying BAM and index..."
ls -lh "${BAM}" "${BAM}.bai"

# ----------------------------------------------------------------------------
# 2. Re-sort BAM by coordinate
# ----------------------------------------------------------------------------
BAM_SORTED="${OUT_DIR}/EP1_chip_coord_sorted.bam"

if [ ! -f "${BAM_SORTED}" ]; then
    echo "[$(date)] Re-sorting BAM by coordinate..."
    samtools sort \
        -@ "${SLURM_CPUS_PER_TASK}" \
        -m 8G \
        -T "${TMP_DIR}/sort_tmp" \
        -o "${BAM_SORTED}" \
        "${BAM}"
    samtools index "${BAM_SORTED}"
    echo "  Done: BAM re-sorted"
else
    echo "  Skipping re-sort: ${BAM_SORTED} already exists"
fi

# ----------------------------------------------------------------------------
# 3. Total mapped reads for CPM normalization
# ----------------------------------------------------------------------------
echo "[$(date)] Counting total mapped reads..."
TOTAL_READS=$(samtools flagstat "${BAM}" | grep "mapped (" | head -1 | awk '{print $1}')
echo "Total mapped reads: ${TOTAL_READS}"
echo "${TOTAL_READS}" > "${OUT_DIR}/total_mapped_reads.txt"

# ----------------------------------------------------------------------------
# 4. Build labeled anchor BED (chr | start | end | loop_class)
#    Sorted to match BAM chromosome order
# ----------------------------------------------------------------------------
echo "[$(date)] Building combined labeled anchor BED..."

declare -A ANCHOR_FILES=(
    ["REH_specific"]="REH_specific_only_05_10kb_exp30kb_anchors_sorted.bed"
    ["EP1_specific"]="EP1_specific_only_05_10kb_exp30kb_anchors_sorted.bed"
    ["ARO_specific"]="ARO_specific_only_05_10kb_exp30kb_anchors_sorted.bed"
    ["shared_EP1_REH"]="shared_EP1_REH_only_05_10kb_exp30kb_anchors_sorted.bed"
    ["shared_ARO_REH"]="shared_ARO_REH_only_05_10kb_exp30kb_anchors_sorted.bed"
    ["shared_EP1_ARO"]="shared_EP1_ARO_only_05_10kb_exp30kb_anchors_sorted.bed"
    ["shared_all3"]="shared_loops_all3_05_10kb_exp30kb_anchors_sorted.bed"
)

LABELED_BED="${OUT_DIR}/all_anchors_labeled.bed"
> "${LABELED_BED}"

for CLASS in "${!ANCHOR_FILES[@]}"; do
    FILE="${ANCHORS_DIR}/${ANCHOR_FILES[$CLASS]}"
    if [ -f "${FILE}" ]; then
        awk -v cls="${CLASS}" 'BEGIN{OFS="\t"} {print $1,$2,$3,cls}' "${FILE}" >> "${LABELED_BED}"
        echo "  Added: ${CLASS} ($(wc -l < "${FILE}") anchors)"
    else
        echo "  WARNING: ${FILE} not found, skipping"
    fi
done

# Sort BED using BAM chromosome order
CHR_ORDER="${TMP_DIR}/chr_order.txt"
samtools view -H "${BAM_SORTED}" \
    | grep "^@SQ" \
    | awk -F'\t' '{print $2}' \
    | sed 's/SN://' \
    > "${CHR_ORDER}"

LABELED_BED_SORTED="${OUT_DIR}/all_anchors_labeled_sorted.bed"

awk 'NR==FNR{rank[$1]=NR; next} {print rank[$1], $0}' \
    "${CHR_ORDER}" "${LABELED_BED}" \
    | sort -k1,1n -k3,3n \
    | awk '{$1=""; sub(/^ /,""); print}' \
    | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' \
    > "${LABELED_BED_SORTED}"

echo "Total anchors in sorted BED: $(wc -l < "${LABELED_BED_SORTED}")"

# ----------------------------------------------------------------------------
# 5. RUN 1 — bedtools coverage default (full stats per anchor)
#
#    Output columns:
#    $1  chr
#    $2  start
#    $3  end
#    $4  loop_class
#    $5  read_count        ← overlapping reads
#    $6  bases_covered     ← number of bases with ≥1 read
#    $7  anchor_length     ← total anchor size (bp)
#    $8  fraction_covered  ← bases_covered / anchor_length
# ----------------------------------------------------------------------------
echo "[$(date)] Running bedtools coverage (full stats)..."

bedtools coverage \
    -a "${LABELED_BED_SORTED}" \
    -b "${BAM_SORTED}" \
    -sorted \
    > "${OUT_DIR}/anchor_coverage_full.bed"

echo "  Done: anchor_coverage_full.bed ($(wc -l < "${OUT_DIR}/anchor_coverage_full.bed") lines)"

# ----------------------------------------------------------------------------
# 6. RUN 2 — bedtools coverage -counts only
#    Faster, lighter — just read_count per anchor, no extra columns
# ----------------------------------------------------------------------------
echo "[$(date)] Running bedtools coverage (counts only)..."

bedtools coverage \
    -a "${LABELED_BED_SORTED}" \
    -b "${BAM_SORTED}" \
    -counts \
    -sorted \
    > "${OUT_DIR}/anchor_read_counts.bed"

# Output: chr | start | end | loop_class | read_count
echo "  Done: anchor_read_counts.bed ($(wc -l < "${OUT_DIR}/anchor_read_counts.bed") lines)"

# ----------------------------------------------------------------------------
# 7. Sanity check
# ----------------------------------------------------------------------------
echo "[$(date)] Per-class summary:"
echo -e "loop_class\tn_anchors\ttotal_reads\tmean_read_count" > "${OUT_DIR}/per_class_summary.tsv"

awk 'BEGIN{OFS="\t"} {
    class[$4]++
    reads[$4]+=$5
} END{
    for(c in class) print c, class[c], reads[c], reads[c]/class[c]
}' "${OUT_DIR}/anchor_read_counts.bed" \
    | sort >> "${OUT_DIR}/per_class_summary.tsv"

cat "${OUT_DIR}/per_class_summary.tsv"

TOTAL_OUT=$(awk 'NR>1 {sum+=$2} END{print sum}' "${OUT_DIR}/per_class_summary.tsv")
echo ""
echo "Total anchors in output: ${TOTAL_OUT} (expected: 11223)"

# ----------------------------------------------------------------------------
# Cleanup
# ----------------------------------------------------------------------------
rm -rf "${TMP_DIR}"

echo "[$(date)] All done! Output in: ${OUT_DIR}"
echo ""
echo "Files for R:"
echo "  anchor_coverage_full.bed  → cols: chr, start, end, loop_class, read_count, bases_covered, anchor_length, fraction_covered"
echo "  anchor_read_counts.bed    → cols: chr, start, end, loop_class, read_count"
echo "  total_mapped_reads.txt    → total mapped reads for CPM: ${TOTAL_READS}"

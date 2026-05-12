#!/bin/bash
#SBATCH --job-name=random_regions_coverage
#SBATCH --output=random_regions_coverage_%j.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00

# ============================================================================
# Random genomic regions (30 kb) x ChIP-seq coverage — bedtools
# Control baseline for E2A-PBX1 ChIP loop anchor analysis
# Designed to match loop_chip_coverage.sh exactly:
#   - same BAM
#   - same region size (30 kb, matching ±10 kb expanded anchors)
#   - same N as original anchor input (11,233)
#   - excludes ENCODE blacklist regions
#   - excludes regions overlapping loop anchors (avoids signal inflation)
# ============================================================================

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
ANCHORS_DIR="/ijc/LABS/STIK/RAW/Hi-C/HiC_marta_Oct25/loops/intersected_loops"
BAM="/ijc/LABS/STIK/DATA/CHIP-seq/public_data_september_2023/bowtie2_results/final_data/697HF_E2APBX1-HF_clean.bam"
OUT_DIR="${ANCHORS_DIR}/chip_coverage"
TMP_DIR="${OUT_DIR}/tmp_random"

CHROM_SIZES="/mnt/beegfs/xvilchez/files/hg38.chrom.sizes.txt"
BLACKLIST="/mnt/beegfs/xvilchez/files/black_list"
ANCHORS_BED="${OUT_DIR}/all_anchors_labeled.bed"

# Match exactly what was used in the ChIP coverage script
N_RANDOM=11233
REGION_SIZE=30000

# Oversample 5x to account for blacklist + anchor exclusion
N_OVERSAMPLE=$(( N_RANDOM * 5 ))

mkdir -p "${OUT_DIR}" "${TMP_DIR}"

# ----------------------------------------------------------------------------
# Modules
# ----------------------------------------------------------------------------
module load BEDTools/2.31.1-GCC-11.2.0
module load SAMtools/1.19.2-foss-2021b

# ----------------------------------------------------------------------------
# 1. Verify inputs
# ----------------------------------------------------------------------------
echo "[$(date)] Verifying input files..."
for f in "${BAM}" "${BAM}.bai" "${CHROM_SIZES}" "${BLACKLIST}" "${ANCHORS_BED}"; do
    if [ ! -f "${f}" ]; then
        echo "  ERROR: Missing file: ${f}"
        exit 1
    fi
done
ls -lh "${BAM}" "${BAM}.bai" "${CHROM_SIZES}" "${BLACKLIST}" "${ANCHORS_BED}"

# ----------------------------------------------------------------------------
# 2. Reuse sorted BAM from ChIP coverage script
# ----------------------------------------------------------------------------
BAM_SORTED="${OUT_DIR}/EP1_chip_coord_sorted.bam"

if [ ! -f "${BAM_SORTED}" ]; then
    echo "[$(date)] Sorted BAM not found — re-sorting..."
    samtools sort \
        -@ "${SLURM_CPUS_PER_TASK}" \
        -m 8G \
        -T "${TMP_DIR}/sort_tmp" \
        -o "${BAM_SORTED}" \
        "${BAM}"
    samtools index "${BAM_SORTED}"
else
    echo "[$(date)] Reusing sorted BAM: ${BAM_SORTED}"
fi

# ----------------------------------------------------------------------------
# 3. Reuse total mapped reads from ChIP coverage script
# ----------------------------------------------------------------------------
READS_FILE="${OUT_DIR}/total_mapped_reads.txt"
if [ -f "${READS_FILE}" ]; then
    TOTAL_READS=$(cat "${READS_FILE}")
    echo "[$(date)] Reusing total mapped reads: ${TOTAL_READS}"
else
    echo "[$(date)] Counting total mapped reads..."
    TOTAL_READS=$(samtools flagstat "${BAM}" | grep "mapped (" | head -1 | awk '{print $1}')
    echo "${TOTAL_READS}" > "${READS_FILE}"
fi

# ----------------------------------------------------------------------------
# 4. Extract chromosome order from BAM header
# ----------------------------------------------------------------------------
echo "[$(date)] Extracting chromosome order from BAM header..."
CHR_ORDER="${TMP_DIR}/chr_order.txt"
samtools view -H "${BAM_SORTED}" \
    | grep "^@SQ" \
    | awk -F'\t' '{print $2}' \
    | sed 's/SN://' \
    > "${CHR_ORDER}"

# Filter chrom.sizes to chromosomes in BAM and standard chr1-22, X, Y only
CHROM_SIZES_FILTERED="${TMP_DIR}/chrom_sizes_filtered.txt"
awk 'NR==FNR{keep[$1]=1; next} keep[$1]' \
    "${CHR_ORDER}" "${CHROM_SIZES}" \
    | grep -E '^chr([0-9]{1,2}|X|Y)\b' \
    > "${CHROM_SIZES_FILTERED}"

echo "  Chromosomes retained: $(wc -l < "${CHROM_SIZES_FILTERED}")"

# ----------------------------------------------------------------------------
# 5. Generate random regions
#    - Oversample 5x (up from 3x) to account for blacklist + anchor exclusion
#    - Fixed seed for reproducibility
#    - Standard chromosomes only (chr1-22, X, Y)
# ----------------------------------------------------------------------------
echo "[$(date)] Generating ${N_OVERSAMPLE} candidate random regions of ${REGION_SIZE} bp..."

RAW_RANDOM="${TMP_DIR}/random_raw.bed"

bedtools random \
    -n "${N_OVERSAMPLE}" \
    -l "${REGION_SIZE}" \
    -g "${CHROM_SIZES_FILTERED}" \
    -seed 42 \
    > "${RAW_RANDOM}"

echo "  Raw random regions generated: $(wc -l < "${RAW_RANDOM}")"

# ----------------------------------------------------------------------------
# 5a. Remove ENCODE blacklist regions
# ----------------------------------------------------------------------------
echo "[$(date)] Removing ENCODE blacklist regions..."

CLEAN_BL="${TMP_DIR}/random_no_blacklist.bed"

bedtools intersect \
    -a "${RAW_RANDOM}" \
    -b "${BLACKLIST}" \
    -v \
    > "${CLEAN_BL}"

N_AFTER_BL=$(wc -l < "${CLEAN_BL}")
echo "  Regions after blacklist exclusion: ${N_AFTER_BL}"

# ----------------------------------------------------------------------------
# 5b. Remove regions overlapping loop anchors
#     Uses the expanded (30 kb) anchor windows so no partial overlaps slip
#     through — this prevents signal inflation in the control set
# ----------------------------------------------------------------------------
echo "[$(date)] Removing regions overlapping loop anchors..."

CLEAN_RANDOM="${TMP_DIR}/random_clean.bed"

bedtools intersect \
    -a "${CLEAN_BL}" \
    -b "${ANCHORS_BED}" \
    -v \
    > "${CLEAN_RANDOM}"

N_AFTER_ANCHORS=$(wc -l < "${CLEAN_RANDOM}")
echo "  Regions after anchor exclusion: ${N_AFTER_ANCHORS} (target: ${N_RANDOM})"

if [ "${N_AFTER_ANCHORS}" -lt "${N_RANDOM}" ]; then
    echo "  ERROR: Not enough clean regions after exclusions."
    echo "         Got ${N_AFTER_ANCHORS}, need ${N_RANDOM}."
    echo "         Increase N_OVERSAMPLE (currently ${N_OVERSAMPLE}) and rerun."
    exit 1
fi

# Trim to exact N_RANDOM
head -n "${N_RANDOM}" "${CLEAN_RANDOM}" > "${TMP_DIR}/random_trimmed.bed"
mv "${TMP_DIR}/random_trimmed.bed" "${CLEAN_RANDOM}"

echo "  Final clean random regions: $(wc -l < "${CLEAN_RANDOM}")"

# ----------------------------------------------------------------------------
# 6. Label and sort to match BAM chromosome order
#    Same sorting logic as loop_chip_coverage.sh
# ----------------------------------------------------------------------------
echo "[$(date)] Labeling and sorting random regions..."

LABELED_RANDOM="${TMP_DIR}/random_labeled.bed"
LABELED_RANDOM_SORTED="${OUT_DIR}/random_regions_labeled_sorted.bed"

awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"random_regions"}' "${CLEAN_RANDOM}" \
    > "${LABELED_RANDOM}"

awk 'NR==FNR{rank[$1]=NR; next} {print rank[$1], $0}' \
    "${CHR_ORDER}" "${LABELED_RANDOM}" \
    | sort -k1,1n -k3,3n \
    | awk '{$1=""; sub(/^ /,""); print}' \
    | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' \
    > "${LABELED_RANDOM_SORTED}"

echo "  Sorted random regions: $(wc -l < "${LABELED_RANDOM_SORTED}") lines"

# ----------------------------------------------------------------------------
# 7. RUN 1 — bedtools coverage full stats
#    Identical command to loop_chip_coverage.sh step 5
#    Output: chr | start | end | loop_class | read_count | bases_covered |
#            anchor_length | fraction_covered
# ----------------------------------------------------------------------------
echo "[$(date)] Running bedtools coverage (full stats)..."

bedtools coverage \
    -a "${LABELED_RANDOM_SORTED}" \
    -b "${BAM_SORTED}" \
    -sorted \
    > "${OUT_DIR}/random_coverage_full_30kb.bed"

echo "  Done: random_coverage_full_30kb.bed ($(wc -l < "${OUT_DIR}/random_coverage_full_30kb.bed") lines)"

# ----------------------------------------------------------------------------
# 8. RUN 2 — bedtools coverage counts only
#    Identical command to loop_chip_coverage.sh step 6
#    Output: chr | start | end | loop_class | read_count
# ----------------------------------------------------------------------------
echo "[$(date)] Running bedtools coverage (counts only)..."

bedtools coverage \
    -a "${LABELED_RANDOM_SORTED}" \
    -b "${BAM_SORTED}" \
    -counts \
    -sorted \
    > "${OUT_DIR}/random_read_counts_30kb.bed"

echo "  Done: random_read_counts_30kb.bed ($(wc -l < "${OUT_DIR}/random_read_counts_30kb.bed") lines)"

# ----------------------------------------------------------------------------
# 9. Sanity check — mirror loop_chip_coverage.sh step 7
# ----------------------------------------------------------------------------
echo "[$(date)] Random regions summary:"
echo -e "loop_class\tn_regions\ttotal_reads\tmean_read_count" \
    > "${OUT_DIR}/random_regions_summary.tsv"

awk 'BEGIN{OFS="\t"} {
    class[$4]++
    reads[$4]+=$5
} END{
    for(c in class) print c, class[c], reads[c], reads[c]/class[c]
}' "${OUT_DIR}/random_read_counts_30kb.bed" \
    | sort >> "${OUT_DIR}/random_regions_summary.tsv"

cat "${OUT_DIR}/random_regions_summary.tsv"

echo ""
echo "[$(date)] Anchor vs random comparison:"
echo "--- Loop anchor classes (from ChIP coverage script) ---"
cat "${OUT_DIR}/per_class_summary.tsv"
echo "--- Random regions ---"
cat "${OUT_DIR}/random_regions_summary.tsv"

# ----------------------------------------------------------------------------
# Exclusion stats summary
# ----------------------------------------------------------------------------
echo ""
echo "[$(date)] Exclusion summary:"
echo "  Raw candidates generated : ${N_OVERSAMPLE}"
echo "  After blacklist removal  : ${N_AFTER_BL}  (removed $(( N_OVERSAMPLE - N_AFTER_BL )))"
echo "  After anchor removal     : ${N_AFTER_ANCHORS}  (removed $(( N_AFTER_BL - N_AFTER_ANCHORS )))"
echo "  Final set used           : ${N_RANDOM}"

# ----------------------------------------------------------------------------
# Cleanup
# ----------------------------------------------------------------------------
rm -rf "${TMP_DIR}"

echo ""
echo "[$(date)] All done! Output in: ${OUT_DIR}"
echo ""
echo "Files for R:"
echo "  random_coverage_full_30kb.bed  → cols: chr, start, end, loop_class, read_count, bases_covered, anchor_length, fraction_covered"
echo "  random_read_counts_30kb.bed    → cols: chr, start, end, loop_class, read_count"
echo "  total_mapped_reads.txt         → shared with anchor script: ${TOTAL_READS}"
echo ""
echo "Note: N_RANDOM=${N_RANDOM} matches the original anchor input."
echo "      Random regions exclude both ENCODE blacklist and all loop anchors"
echo "      (${ANCHORS_BED}) to prevent signal inflation in the control set."

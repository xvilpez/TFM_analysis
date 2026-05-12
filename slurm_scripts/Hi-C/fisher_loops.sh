#!/bin/bash
#SBATCH --job-name=fisher_loops
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=fisher_loops_%A_%a.log

module load BEDTools/2.31.1-GCC-11.2.0

# Paths
WORKDIR="/ijc/LABS/STIK/RAW/Hi-C/HiC_marta_Oct25/loops/intersected_loops"
CHIP=/mnt/beegfs/xvilchez/chip_publicdata/top25_697HF_E2APBX1-HF_peaks.bed
GENOME=/mnt/beegfs/xvilchez/files/hg38.chrom.sizes.txt

cd "$WORKDIR" || exit 1

# Sort ChIP file once
CHIP_SORTED="${WORKDIR}/chip_top25_sorted.bed"
bedtools sort -i "$CHIP" -g "$GENOME" > "$CHIP_SORTED"

for ANCHORS in *_anchors.bed; do
    BASENAME=$(basename "$ANCHORS" _anchors.bed)

    echo "Processing $ANCHORS..."

    # Ensure anchors are sorted & unique (important for fisher)
    bedtools sort -i "$ANCHORS" -g "$GENOME" | uniq > "${BASENAME}_anchors_sorted.bed"

    # Run bedtools fisher
    bedtools fisher \
        -a "${BASENAME}_anchors_sorted.bed" \
        -b "$CHIP_SORTED" \
        -g "$GENOME" > "${BASENAME}_top25_fisher.txt"

    echo "Finished $ANCHORS → ${BASENAME}_top25_fisher.txt"
done

# Cleanup
rm "$CHIP_SORTED"

echo "All jobs completed."


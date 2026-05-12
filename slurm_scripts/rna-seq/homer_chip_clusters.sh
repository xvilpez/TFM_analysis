#!/bin/bash
#SBATCH --job-name=homer_chip_clusters
#SBATCH --output=homer_chip_%A_%a.out
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=1-8

# Load HOMER
module load homer/5.1

# Go to working directory
cd /ijc/LABS/STIK/RAW/rna_seq/rnaseq_marta_REHAro_dec25/DEG_results/Clusters/homer_chip_motifs

# Cluster index from SLURM array
CLUSTER=${SLURM_ARRAY_TASK_ID}

# Input BED file
BEDFILE="cluster_${CLUSTER}_chip_peaks.bed"

# Output directory
OUTDIR="cluster_${CLUSTER}_chip_DEGbg"

# Run HOMER
findMotifsGenome.pl ${BEDFILE} hg38 ${OUTDIR} \
    -size given \
    -mask \
    -bg background_deg_chip_peaks.bed \
    -p ${SLURM_CPUS_PER_TASK}

echo "Finished Cluster ${CLUSTER}"

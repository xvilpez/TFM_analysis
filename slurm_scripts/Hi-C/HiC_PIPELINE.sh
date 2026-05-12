#!/bin/sh
#SBATCH --job-name=HICpipeline
#SBATCH --mem=100gb
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=HICpipeline_%A-%a.log

# ========== VARIABLES ==========
# put in a file call samples.txt the name of the variables

describer=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

source ./config.sh

for dir in "${path_bam}" ; do
  if [ ! -d "${dir}" ]; then
    mkdir -p "${dir}"
  fi
done

# ========== MODULES ==========
module load fastqc-0.11.9-gcc-11.2.0-dd2vd2m
module load Trim_Galore/0.6.6-foss-2021b-Python-3.8.5
module load Perl/5.34.0-GCCcore-11.2.0
module load  R/4.2.1-foss-2021b
module load Bowtie2/2.4.4.1-GCC-11.2.0

# ========== FASTQC ==========
echo "................................................................ 1. START_FASTQC ${describer} ................................................................"

#fastqc ${path_fq}/${describer}_*.fastq.gz -o ${path_fq}

echo "................................................................ 1. END_FASTQC ${describer} ................................................................"

# ========== TRIMMING ==========
echo "................................................................ 2. START_TRIM_GALORE ${describer} ................................................................"

#trim_galore --output_dir ${path_fq}  --paired ${path_fq}/${describer}_*R1*.fastq.gz ${path_fq}/${describer}_*R2*.fastq.gz

echo "................................................................ 2. END_TRIM_GALORE ${describer} ................................................................"

echo "................................................................ 3. START_HICUP_TRUNCATED ${describer} ................................................................"

#perl ${HICUP_trunc} --re1 ${restriction_enzyme} ${path_fq}/${describer}_*R1*_val_1.fq.gz ${path_fq}/${describer}_*R2*_val_2.fq.gz --zip --outdir ${path_fq}

echo "................................................................ 3. END_HICUP_TRUNCATED ${describer} ................................................................"


# ========== ALIGNMENT ==========
echo "................................................................ 4. START_R1_BOWTIE2 ${describer} ................................................................"
bowtie2 --local -x ${indexgenome} --threads 8 \
    -U ${path_fq}/${describer}_*R1*_val_*.trunc.* --reorder -S ${path_bam}/${describer}_R1.sam

echo "................................................................ 4. END_R1_BOWTIE2 ${describer} ................................................................"

echo "................................................................ 5. START_R2_BOWTIE2 ${describer} ................................................................"

bowtie2 --local -x ${indexgenome} --threads 8 \
    -U ${path_fq}/${describer}_*R2*_val_*.trunc.* --reorder -S ${path_bam}/${describer}_R2.sam

echo "................................................................ 5. END_R2_BOWTIE2 ${describer} ................................................................"

# ==========  LAUNCH ANALYSIS SCRIPTS ==========
# Only submit once from task 1 (or use a separate submission script)
if [ "${SLURM_ARRAY_TASK_ID}" == "1" ]; then
    # Get the current job ID (the alignment job array)
    ALIGNMENT_JOB=${SLURM_ARRAY_JOB_ID}
    
    echo "Submitting downstream analysis with dependencies..."
    
    # Submit tagdir to run AFTER all alignment tasks complete
    sbatch --dependency=afterok:${ALIGNMENT_JOB} --array=1-${N} scripts/tagdir.sh
    echo "Tagdir jobs submitted (depends on ${ALIGNMENT_JOB})"
    
    # Submit hicExplorer to run AFTER all alignment tasks complete (parallel with tagdir)
    sbatch --dependency=afterok:${ALIGNMENT_JOB} --array=1-${N} scripts/hicExplorer_analysis.sh
    echo "HicExplorer jobs submitted (depends on ${ALIGNMENT_JOB})"
fi

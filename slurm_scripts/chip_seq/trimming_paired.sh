#!/bin/sh
#SBATCH --job-name=trimming
#SBATCH --mem=70gb
#SBATCH --time=11:00:00
#SBATCH --output=trimming_%A-%a.log
#SBATCH --array=1-3

module load Trim_Galore/0.6.6-foss-2021b-Python-3.8.5

names=( REH_CTCF_rep1
REH_CTCF_rep2
REH_INPUT 
) 
describer=${names[$SLURM_ARRAY_TASK_ID-1]}

WD=/ijc/LABS/STIK/RAW/chip_seq/Xavi_TFM/REH_public_marta

mkdir -p ${WD}/Trim_Galore

cd ${WD}/fastq_files

trim_galore --fastqc --output_dir ../Trim_Galore --paired ${describer}_1.fastq.gz ${describer}_2.fastq.gz

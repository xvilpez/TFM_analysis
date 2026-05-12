#!/bin/sh
#SBATCH --job-name=trimming
#SBATCH --mem=70gb
#SBATCH --time=11:00:00
#SBATCH --output=trimming_%A-%a.log
#SBATCH --array=1-1

module load Trim_Galore/0.6.6-foss-2021b-Python-3.8.5

names=( CD19.RO_01736 ) 
describer=${names[$SLURM_ARRAY_TASK_ID-1]}

mkdir -p ../data/Trim_Galore

trim_galore --fastqc --output_dir ../data/Trim_Galore ../data/raw/${describer}.fastq.gz 



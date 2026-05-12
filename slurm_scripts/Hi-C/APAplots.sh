#!/bin/bash

#SBATCH --job-name=aggregate_matrix
#SBATCH --mem=180gb
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=agg_all_%A-%a.log
#SBATCH --array=1-7

module load HiCExplorer/3.7.2-foss-2021b
module load cooler/0.9.1-foss-2021b
module load scikit-learn/0.24.2-foss-2021b
#module load python/3.8.5

names=( REH HiC-REH-EP1 HiC-REH-EP1-1 HiC-REH-EP1-2 HiC-REH-EP1-Aro HiC-REH-EP1-Aro-1 HiC-REH-EP1-Aro-2 )
describer=${names[$SLURM_ARRAY_TASK_ID-1]}

path_mat='/ijc/LABS/STIK/RAW/Hi-C/HiC_marta_Oct25/coolMatrix'
path_out='/ijc/LABS/STIK/RAW/Hi-C/HiC_marta_Oct25/APAplots'
path_bed='/mnt/beegfs/xvilchez/chip_publicdata'

mkdir -p ${path_out}

for peaks in 697HF_E2APBX1-HF_peaks
do

echo "------------------------------------------start ${describer} ${peaks} -----------------------------------------------"

hicAggregateContacts -m ${path_mat}/${describer}_10kb_KR.cool \
     --outFileName ${path_out}/scaled/${describer}_${peaks}_APA.pdf \
     --BED ${path_bed}/${peaks}.bed \
     --numberOfBins 51 --vMin 1 --vMax 1.3 \
     --largeRegionsOperation center --range 500000:10000000 \
     --transform obs/exp --outFilePrefixMatrix ${path_out}/scaled/${describer}_${peaks} \
     --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX --operationType mean --mode intra-chr

echo "------------------------------------------end ${describer} ${peaks} -----------------------------------------------"

done


#!/bin/sh
#SBATCH --job-name=all_chip_pipeline
#SBATCH --mem=80gb
#SBATCH --time=20:00:00
#SBATCH --output=chip_pip_%A-%a.log
#SBATCH --cpus-per-task=8
#SBATCH --array=1-3

module load Bowtie2/2.4.4.1-GCC-11.2.0
module load SAMtools/1.13-foss-2021b
module load picard/2.26.3-Java-11
module load bedtools2-2.30.0-gcc-11.2.0-q7z4zez
module load deepTools/3.5.1-foss-2021b

names=( REH_CTCF_rep1
REH_CTCF_rep2
REH_INPUT
)

describer=${names[$SLURM_ARRAY_TASK_ID-1]}

WD=/ijc/LABS/STIK/RAW/chip_seq/Xavi_TFM/REH_public_marta

mkdir -p ${WD}/bowtie2_results
mkdir -p ${WD}/bowtie2_results/temp_data
mkdir -p ${WD}/bowtie2_results/final_data
mkdir -p ${WD}/bigwig

#echo "############################################################ 0. alignment ${describer} ###########################################################################"
#bowtie2 -q --very-sensitive-local -x /mnt/beegfs/public/references/index/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as \
#    -1 ${WD}/Trim_Galore/${describer}_1_val_1.fq.gz -2 ${WD}/Trim_Galore/${describer}_2_val_2.fq.gz\
#    --threads 8 \
#    --local    \
#    -S ${WD}/bowtie2_results/${describer}.sam

#echo "############################################################ 0. alignment ${describer} ###########################################################################"


echo "############################################################ 1. START sam -> bam + sort + index ${describer} ###########################################################################"

samtools view -bS ${WD}/bowtie2_results/${describer}.sam | samtools sort -O bam -o ${WD}/bowtie2_results/${describer}.bam -T ${WD}/bowtie2_results/tmp_${describer}


echo "############################################################# 1. END sam -> bam + sort + index ${describer} ###########################################################################"

echo "############################################################ 2. START_remove_mtDNA_READS $describer ###################################################################################"

samtools view -h ${WD}/bowtie2_results/${describer}.bam | grep -v chrM | samtools sort -O bam -o ${WD}/bowtie2_results/temp_data/${describer}.rmChrM.bam -T ${WD}/bowtie2_results/temp_data

echo "########################################################### 2. END_remove_mtDNA_READS $describer #######################################################################################"

echo "############################################################# 3. START_filtering_protperly_reads $describer ###########################################################################"

samtools view -F 2304 -b -q 10  ${WD}/bowtie2_results/temp_data/${describer}.rmChrM.bam | samtools sort -O bam -o ${WD}/bowtie2_results/temp_data/${describer}.qual.bam -T ${WD}/bowtie2_results/temp_data

# -b: sortida en bam
# -q 10: filtra les lectures si es menor a 10 fora

echo "########################################################### 3. END_filtering_protperly_reads $describer ###############################################################################"

echo "############################################################## 4. START_mark_duplicates $describer #################################################################"

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
I=${WD}/bowtie2_results/temp_data/${describer}.qual.bam \
O=${WD}/bowtie2_results/temp_data/${describer}_removed_duplicates.bam \
M=${WD}/bowtie2_results/temp_data/${describer}_marked_dup_metrics.txt \
REMOVE_DUPLICATES=true ASSUME_SORTED=true VERBOSITY=WARNING


echo "############################################################ 4. END_mark_duplicates $describer #####################################################################"

echo "############################################################ 5. START_blacklist regions + index $describer ############################################################################"

bedtools intersect -nonamecheck -v -abam ${WD}/bowtie2_results/temp_data/${describer}_removed_duplicates.bam -b /mnt/beegfs/xvilchez/files/black_list > ${WD}/bowtie2_results/final_data/${describer}_clean.bam
samtools index ${WD}/bowtie2_results/final_data/${describer}_clean.bam

echo "############################################################ 5. END_blacklist + index $describer ######################################################################################"



echo "############################################################## 8. START_bamcoverage $describer ####################################################################"

bamCoverage --bam ${WD}/bowtie2_results/final_data/${describer}_clean.bam \
        --outFileName ${WD}/bigwig/${describer}.bw \
        --effectiveGenomeSize 2913022398 \
        --outFileFormat bigwig \
        --binSize 1 --normalizeUsing RPGC > ${WD}/bigwig/${describer}.log

echo "############################################################## 8. END_bamcoverage $describer ##########################################################################"

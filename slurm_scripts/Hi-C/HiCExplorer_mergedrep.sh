#!/bin/sh
#SBATCH --job-name=hicMatrix_mergedreplicates
#SBATCH --mem=150gb
#SBATCH --time=30:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=hicMatrix_mergedrep%A-%a.log
#SBATCH --array=1-1

# ========== VARIABLES ==========
# put in a file call samples.txt the name of the variables

source ./config.sh

names=( HiC-REH-EP1-Aro )

describer=${names[$SLURM_ARRAY_TASK_ID-1]}

module load HiCExplorer/3.7.6-foss-2021b

echo "................................................................ START hicSumMatrices and hicCorrectMatrix and hicConvertFormat for ${describer} ................................................................"
for number in 5 10 20 50 100
do
    hicSumMatrices --matrices ${path_hicMatrix}/${describer}-1_${number}kb.h5 ${path_hicMatrix}/${describer}-2_${number}kb.h5 --outFileName ${path_hicMatrix}/${describer}_${number}kb.h5

    echo "................................................................ START hicCorrect KR ${number}kb ${describer} ................................................................"
    
    hicCorrectMatrix correct --correctionMethod KR \
         --matrix ${path_hicMatrix}/${describer}_${number}kb.h5  \
         --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
         --perchr --outFileName ${path_hicMatrix}/${describer}_${number}kb_KR.h5

    echo "................................................................ END hicCorrect KR ${number}kb ${describer} ................................................................"

    echo "................................................................ START hicConvert to cool ${number}kb ${describer} ................................................................"
    hicConvertFormat --matrices ${path_hicMatrix}/${describer}_${number}kb_KR.h5 \
        --outFileName ${path_coolMatrix}/${describer}_${number}kb_KR.cool \
        --inputFormat h5 \
        --outputFormat cool
    echo "................................................................ END hicConvert to cool ${number}kb ${describer} ................................................................"

done

# ==========  LOOP CALLING ==========
module purge
unset PYTHONPATH
unset LD_LIBRARY_PATH

module load Miniconda3/24.7.1-0
source /opt/ohpc/pub/libs/easybuild/4.5.0/software/Miniconda3/24.7.1-0/etc/profile.d/conda.sh
conda activate mustache

for number in 5 10 20
do

    echo " ................................................................ START mustache ${describer} ................................................................ "
    python -m mustache -f ${path_coolMatrix}/${describer}_${number}kb_KR.cool  -r ${number}kb -pt 0.1 -o ${path_loops}/hic_${describer}_loops_01_${number}kb.tsv
    python -m mustache -f ${path_coolMatrix}/${describer}_${number}kb_KR.cool  -r ${number}kb -pt 0.05 -o ${path_loops}/hic_${describer}_loops_05_${number}kb.tsv
    echo " ................................................................ END mustache ${describer} ................................................................ "

done

conda deactivate

module purge
module load cooltools/0.5.2-foss-2021b
module load cooler



# ==========  GENERATE SADDLE PLOTS ==========
echo " ................................................................ START cooltools expected-cis ${describer} ................................................................ "
cooltools expected-cis -p 8 -o ${path_cooltools}/${describer}_100kb_KR_exp.tsv \
   ${path_coolMatrix}/${describer}_100kb_KR.cool
echo " ................................................................ END cooltools expected-cis ${describer} ................................................................ "



echo " ................................................................ START cooltools eigs-cis ${describer} ................................................................ "
cooltools eigs-cis ${path_coolMatrix}/${describer}_100kb_KR.cool \
--phasing-track ${ref_compartments} \
     -o ${path_cooltools}/${describer}_100kb_KR_ev_ac
echo " ................................................................ END cooltools eigs-cis ${describer} ................................................................ "

echo " ................................................................ START awk ${describer} ................................................................ "
awk 'NR>1 {print $1, $2, $3, $5}' ${path_cooltools}/${describer}_100kb_KR_ev_ac.cis.vecs.tsv > ${path_cooltools}/${describer}_100kb_EV1_ac.bedgraph
awk 'NR>1 {print $1, $2, $3, $6}' ${path_cooltools}/${describer}_100kb_KR_ev_ac.cis.vecs.tsv > ${path_cooltools}/${describer}_100kb_EV2_ac.bedgraph
awk 'NR>1 {print $1, $2, $3, $7}' ${path_cooltools}/${describer}_100kb_KR_ev_ac.cis.vecs.tsv > ${path_cooltools}/${describer}_100kb_EV3_ac.bedgraph
echo " ................................................................ END awk ${describer} ................................................................ "

echo " ................................................................ START cooltools saddle ${describer} ................................................................ "

cooltools saddle --qrange 0.02 0.98 --strength \
    --vmin 0.2 --vmax 4 \
    -o ${path_cooltools}/${describer} --fig pdf \
   ${path_coolMatrix}/${describer}_100kb_KR.cool \
    ${path_cooltools}/${describer}_100kb_KR_ev_ac.cis.vecs.tsv \
    ${path_cooltools}/${describer}_100kb_KR_exp.tsv

echo " ................................................................ END cooltools saddle ${describer} ................................................................ "

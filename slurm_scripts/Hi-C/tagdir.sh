#!/bin/sh

#SBATCH --job-name=homer
#SBATCH --mem=60gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=tag_dir_%A-%a.txt

# ========== VARIABLES ==========
# put in a file call samples.txt the name of the variables

describer=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)
source ./config.sh

for dir in "${path_homer}" ; do
  if [ ! -d "${dir}" ]; then
    mkdir -p "${dir}"
  fi
done

# ========== MODULES ==========

module load homer/0.1
module load Java/17.0.2



# echo " ................................................................ START makeTagDirectory 1 ${describer} ................................................................"

# makeTagDirectory ${path_homer}/${describer}_unfiltered \
#                  ${path_bam}/${describer}_R1.sam,${path_bam}/${describer}_R2.sam \
#                  -tbp 1 -illuminaPE

# echo " ................................................................ END makeTagDirectory 1 ${describer} ................................................................"

echo " ................................................................ START cp ${describer} ................................................................ "
cp -r ${path_homer}/${describer}_unfiltered ${path_homer}/${describer}_filtered
echo " ................................................................ END cp ${describer} ................................................................"

echo " ................................................................ START makeTagDirectory 2 ${describer} ................................................................ "

makeTagDirectory ${path_homer}/${describer}_filtered -update \
                -genome ${genome} -removePEbg \
                -restrictionSite ${restrictionSequence} -both \
                -removeSelfLigation -removeSpikes 10000 5

echo " ................................................................ END makeTagDirectory 2 ${describer} ................................................................ "


echo "................................................................ START hic file ${describer} ................................................................"

tagDir2hicFile.pl ${path_homer}/${describer}_filtered  -juicer auto -genome ${genome} -juicerExe "java -jar /ijc/LABS/STIK/RAW/xavi/juicer_tools.1.9.9_jcuda.0.8.jar" -p 8

echo "................................................................ END hic file ${describer} ................................................................"


# Only task 1 submits the next job
if [ "${SLURM_ARRAY_TASK_ID}" == "1" ]; then
    # Get the current tagdir job ID
    TAGDIR_JOB=${SLURM_ARRAY_JOB_ID}

    echo "Submitting txtfile analysis with dependency..."

    # Submit txtfile to run AFTER all tagdir tasks complete
    sbatch --dependency=afterok:${TAGDIR_JOB} --array=1-${N} scripts/txtfile.sh
    echo "Txtfile jobs submitted (depends on ${TAGDIR_JOB})"
fi

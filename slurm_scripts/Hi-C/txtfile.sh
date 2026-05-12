#!/bin/bash

#SBATCH --job-name=txt_homer
#SBATCH --mem=60gb
#SBATCH --time=26:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=homertxt_%a.txt


# ========== VARIABLES ==========
# put in a file call samples.txt the name of the variables

describer=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)
source ./config.sh

# ========== MODULES ==========

module load homer/0.1
module load Java/17.0.2

echo " ................................................................ START .txt file ${describer} ................................................................ "

tagDir2hicFile.pl ${path_homer}/${describer}_filtered -genome ${genome} -juicerExe "java -jar /ijc/LABS/STIK/RAW/xavi/juicer_tools.1.9.9_jcuda.0.8.jar" -p 8

echo " ................................................................ END .txt file ${describer} ................................................................ "

# Only task 1 submits the next job
if [ "${SLURM_ARRAY_TASK_ID}" == "1" ]; then
    # Get the current txtfile job ID
    TXTFILE_JOB=${SLURM_ARRAY_JOB_ID}

    echo "Submitting cscore analysis with dependency..."

    # Submit cscore to run AFTER all txtfile tasks complete
   #sbatch --dependency=afterok:${TXTFILE_JOB} --array=1-${N} scripts/cscore.sh
    sbatch --dependency=afterok:${TXTFILE_JOB} --array=1-${N} scripts/cscore_m.sh
    echo "Cscore jobs submitted (depends on ${TXTFILE_JOB})"
fi

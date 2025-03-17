#!/bin/bash

# Definining resource we want to allocate.
#SBATCH --nodes=1
#SBATCH --ntasks=40

# 5 CPU cores per task to keep the parallel data feeding going. 
#SBATCH --cpus-per-task=1

# Allocate enough memory.
#SBATCH --mem=48G
#SBATCH -p small
###SBATCH -p gputest

# Time limit on Puhti's small partition is 3 days. 72:00:00
#SBATCH -t 72:00:00
###SBATCH -t 00:30:00
#SBATCH -J proc_trig

# Puhti project number
#SBATCH --account=Project_2001426

# Log file locations, %j corresponds to slurm job id. symlinks didn't work. Will add hard links to directory instead. Now it saves in projappl dir.
#SBATCH -o logs/%j.out
#SBATCH -e logs/%j.err

FOLDER_ID="$1"
TRIGGER_INPUT_DIR="$2"
RE_OUTPUT_DIR="$3" 
SAME_ORG_EID_DIR="$4" 
JOB_DIR="$5"

# delete job file on exit
function on_exit {
    rm -f $JOB_DIR/$SLURM_JOBID
}
trap on_exit EXIT

source /scratch/project_2001426/qingyao/new-trigger-out-test/.pyvenv/bin/activate

for f in ${TRIGGER_INPUT_DIR}/${FOLDER_ID}/*.tar.gz; do
    subfolder_id=$(basename $f .tar.gz)
    python3 filter_relationship.py \
            --in_tar_gz_ann_path $f \
            --in_RE_outtsvgz_path ${RE_OUTPUT_DIR}/output_${FOLDER_ID}/${subfolder_id}.outtsv.gz  \
            --out_filtered_relation_with_eid ${SAME_ORG_EID_DIR}/${FOLDER_ID}/same-org-outputs-with-eids.tsv.gz
done

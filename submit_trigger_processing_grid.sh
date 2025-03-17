#!/bin/bash
MAX_JOBS=40

FOLDER_ID_START="$1" 
FOLDER_ID_END="$2" 
TRIGGER_INPUT_DIR="$3" 
RE_OUTPUT_DIR="$4" 
SAME_ORG_EID_DIR="$5" 
JOB_DIR="$6"

for i in $(seq ${FOLDER_ID_START} ${FOLDER_ID_END}); do
    while true; do
        jobs=$(ls "$JOB_DIR" | wc -l)
        if [ $jobs -lt $MAX_JOBS ]; then break; fi
        echo "Too many jobs, sleeping ..."
        sleep 60
    done
    ## about 4 min per job
    job_id=$(
        sbatch submit_trigger_processing_single.sh ${i} ${TRIGGER_INPUT_DIR} ${RE_OUTPUT_DIR} ${SAME_ORG_EID_DIR} ${JOB_DIR}\
        | perl -pe 's/Submitted batch job //'
    )
    touch ${JOB_DIR}/$job_id
    echo "Submitted batch job $job_id"
    sleep 10
done
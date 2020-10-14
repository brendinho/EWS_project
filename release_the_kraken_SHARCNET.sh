#!/bin/bash

module load gcc/7.3.0

mkdir -p Jobs
rm Skynet_hesitance

GRAHAM_MAX=1000

MAX_HOURS=23
MAX_MINS=59

PARAMETER_FILE_NAME="job_table.txt"

TABLE_LENGTH=`cat $PARAMETER_FILE_NAME | wc -l`
if [ $TABLE_LENGTH -eq 0 ]; then exit 1; fi

NUM_JOBS_QUEUED=`squeue -u b2philli | wc -l`
QUEUE_SLACK=`expr $GRAHAM_MAX - $NUM_JOBS_QUEUED - 1`
if [ $QUEUE_SLACK -eq 0 ]; then exit 1; fi

make

MAX_TRIALS="$(($TABLE_LENGTH>$QUEUE_SLACK ? $QUEUE_SLACK : $TABLE_LENGTH))"

for i in `seq 1 $MAX_TRIALS`; do
	sbatch --time=0-$MAX_HOURS:$MAX_MINS individual_job_SHARCNET.sh $PARAMETER_FILE_NAME $i $MAX_HOURS $(($MAX_MINS-10))
	sleep 1.0
done

sed -i -e 1,$(($MAX_TRIALS))d $PARAMETER_FILE_NAME

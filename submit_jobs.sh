# !/bin/bash

mkdir -p Jobs
rm Skynet_hesitance

make

MAX_HOURS=23
MAX_MINS=59

PARAMETER_FILE_NAME="job_table.txt"

HAGRID_MAX=95
TABLE_LENGTH=1

NUM_JOBS_RUNNING_HAGRID=`squeue -u b2philli | grep hagrid | wc -l` # -t RUNNING
QUEUE_SLACK_HAGRID=`expr $HAGRID_MAX - $NUM_JOBS_RUNNING_HAGRID`

while [ "$QUEUE_SLACK_HAGRID" -ge 1 ]; do

	TABLE_LENGTH=`cat job_table.txt | wc -l`

	if [ "$TABLE_LENGTH" -eq 0 ]; then
        break
    fi

	MAX_TRIALS="$(($TABLE_LENGTH>$QUEUE_SLACK_HAGRID ? $QUEUE_SLACK_HAGRID : $TABLE_LENGTH))"

	for i in `seq 1 $MAX_TRIALS`; do
		sbatch --time=0-$MAX_HOURS:$MAX_MINS individual_job_HAGRID.sh $PARAMETER_FILE_NAME $i $MAX_HOURS $(($MAX_MINS-10))
		sleep 0.5
	done

	sleep 10

	sed -i -e 1,$(($MAX_TRIALS))d job_table.txt

	NUM_JOBS_RUNNING_HAGRID=`squeue -u b2philli | grep hagrid | wc -l` # -t RUNNING
	QUEUE_SLACK_HAGRID=`expr $HAGRID_MAX - $NUM_JOBS_RUNNING_HAGRID`

done


NORMAL_MAX=60
TABLE_LENGTH=1

NUM_JOBS_RUNNING_NORMAL=`squeue -u b2philli | grep hpc | wc -l` # -t RUNNING
QUEUE_SLACK_NORMAL=`expr $NORMAL_MAX - $NUM_JOBS_RUNNING_NORMAL`

while [ "$QUEUE_SLACK_NORMAL" -ge 1 ]; do

	TABLE_LENGTH=`cat job_table.txt | wc -l`

	if [ "$TABLE_LENGTH" -eq 0 ]; then
        break
    fi

	MAX_TRIALS="$(($TABLE_LENGTH>$QUEUE_SLACK_NORMAL ? $QUEUE_SLACK_NORMAL : $TABLE_LENGTH))"

	for i in `seq 1 $MAX_TRIALS`; do
		sbatch --time=0-$MAX_HOURS:$MAX_MINS individual_job_NORMAL.sh $PARAMETER_FILE_NAME $i $MAX_HOURS $(($MAX_MINS-10))
		sleep 0.5
	done

	sleep 10

	sed -i -e 1,$(($MAX_TRIALS))d job_table.txt

	NUM_JOBS_RUNNING_NORMAL=`squeue -u b2philli | grep hpc | wc -l` # -t RUNNING
	QUEUE_SLACK_NORMAL=`expr $NORMAL_MAX - $NUM_JOBS_RUNNING_NORMAL`

done

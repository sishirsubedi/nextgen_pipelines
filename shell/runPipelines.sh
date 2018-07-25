#!/usr/bin/env bash
SHELL=/bin/bash

currentdate=`date '+%Y-%m-%d %H:%M:%S'`
echo "Running pipelines cron job at - $currentdate "

###database credentials --
user=hhadmin
password=ngs3127
database=test

################################################################################
# running pipeline from sampleAnalysisQueue
################################################################################

statement="select queueID from sampleAnalysisQueue where status=0 order by queueID asc limit 5;"
jobs=()

while  read -r queueID ;
do
   echo "queueID is $queueID"
   jobs+=($queueID)
done < <(mysql --user="$user" --password="$password" --database="$database" --execute="$statement" -N)



## run jobs in parallel
for job in "${jobs[@]}"
do
    echo "running pipeline for queueID- $job"
    nohup bash /home/pipelines/master/shell/pipelineThread.sh -q $job  &
    #bash /home/pipelines/master/shell/pipelineThread.sh -q $job
    pids[${i}]=$!
done


# ##### wait for all background process to end
# for pid in ${pids[*]}; do
#     echo "pid for queue $queueID is  $pid"
#     wait $pid
# done

echo "cron job finished running at - $currentdate"

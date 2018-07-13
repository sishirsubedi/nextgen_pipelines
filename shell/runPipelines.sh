#!/usr/bin/env bash
SHELL=/bin/bash

###database credentials
user=root
password=molSeq3127
database=test
table=sampleAnalysisQueue

################################################################################
# running pipeline from sampleAnalysisQueue
################################################################################

defaultstatus=0
statement="select queueID from $table where status=$defaultstatus order by queueID asc limit 5;"
jobs=()

while  read -r queueID ;
do
   echo "queueID is $queueID"
   jobs+=($queueID)
done < <(mysql --user="$user" --password="$password" --database="$database" --execute="$statement" -N)

## run jobs in parallel
for job in "${jobs[@]}"
do
    echo "job is $job"
    nohup bash /home/pipelines/master/shell/pipelineThread.sh -q $job  &

    pids[${i}]=$!

done


##### wait for all background process to end
for pid in ${pids[*]}; do
    echo "pid is  $pid"
    wait $pid
done


################################################################################
# updating analysis results
################################################################################

# update analysis
analysischeck_statement="select queueID,plStatus from pipelineStatus where plStatus='runCompleted' order by queueID asc limit 5;"
while  read -r queueID plStatus;
do
  newStatus='addingAnalysis'
  echo "Running analysis"
  updateanalysis_statement="update pipelineStatus set plStatus='$newStatus' where queueID=$queueID and plStatus='$plStatus'"
  mysql --user="$user" --password="$password" --database="$database" --execute="$updateanalysis_statement"
  bash /home/pipelines/master/shell/runAnalysis.sh -q $queueID
done< <(mysql --user="$user" --password="$password" --database="$database" --execute="$analysischeck_statement" -N)


DATE=`date '+%Y-%m-%d %H:%M:%S'`
echo "cron job finished running at - $DATE" >> /home/cron_log.txt

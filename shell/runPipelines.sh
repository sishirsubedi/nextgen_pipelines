#!/usr/bin/env bash
SHELL=/bin/bash

###database credentials
user=root
password=molSeq3127
database=test


################################################################################
# running pipeline from sampleAnalysisQueue
################################################################################

table=sampleAnalysisQueue
defaultstatus=0
## currently first come first serve-- one at a time
statement="select queueID,  runID,  sampleID, coverageID, vcallerID,  assayID, instrumentID, environmentID, status from $table where status=$defaultstatus order by queueID asc limit 1;"
# IFS='@'
while  read -r queueID  runID  sampleID coverageID vcallerID assayID instrumentID environmentID status;
do

  # first update sampleAnalysisQueue table and set status of this queue to 1 i.e started processing
  updatestatement="UPDATE $table SET status=1 WHERE queueID = $queueID;"
  mysql --user="$user" --password="$password" --database="$database" --execute="$updatestatement"

  # second insert into pipelineStatus table to update status as started
  insertstatement="INSERT INTO pipelineStatus (queueID, plStatus, timeUpdated) VALUES ('$queueID','started',now());"
  mysql --user="$user" --password="$password" --database="$database" --execute="$insertstatement"

  # run pipeline

  echo " $queueID:$runID:$sampleID:$status"

  if [ "$instrumentID" == "proton" ] || [ "$instrumentID" == "pgm" ]
  then
     echo "proton or pgm"
     bash /home/pipelines/master/shell/ionPipelineInterface.sh -r $runID -s $sampleID -c $coverageID -v $vcallerID -a $assayID -i $instrumentID -e $environmentID -q $queueID
  elif [ "$instrumentID" == "nextseq" ] || [ "$instrumentID" == "miseq" ]
  then
     echo "nextseq or miseq"
     bash /home/pipelines/master/shell/illuminaPipelineInterface.sh -r $runID -s $sampleID -a $assayID -i $instrumentID -e $environmentID #-q $queueID
  fi

done< <(mysql --user="$user" --password="$password" --database="$database" --execute="$statement" -N)


################################################################################
# updating analysis results
################################################################################

# update analysis
analysischeck_statement="select queueID,plStatus from pipelineStatus where plStatus='runCompleted' order by queueID asc limit 1;"
while  read -r queueID plStatus;
do
  newStatus='addingAnalysis'
  updateanalysis_statement="update pipelineStatus set plStatus='$newStatus' where queueID=$queueID and plStatus='$plStatus'"
  mysql --user="$user" --password="$password" --database="$database" --execute="$updateanalysis_statement"
  bash /home/pipelines/master/shell/runAnalysis.sh -q $queueID
done < <(mysql --user="$user" --password="$password" --database="$database" --execute="$analysischeck_statement" -N)


DATE=`date '+%Y-%m-%d %H:%M:%S'`
echo "cron job finished running at - $DATE" >> /home/cron_log.txt

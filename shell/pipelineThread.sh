#!/usr/bin/env bash
SHELL=/bin/bash

if [ $# -eq 0 ]
then
	echo "Usage: pipelineThread.sh"
	echo "-q queueID"
	exit
fi

if test $# -gt 0
	then
	while getopts :q: opt
	do
	case $opt in
  q)
		queueID=$OPTARG
		;;
	:)
		echo "Option -$OPTARG requires an argument."
		;;
	\?)
		echo "Invalid option: -$OPTARG"
	esac
	done
	shift $((OPTIND-1))
fi

###database credentials
user=root
password=molSeq3127
database=test
table=sampleAnalysisQueue


################################################################################
# running pipeline from sampleAnalysisQueue
################################################################################

job_statement="select queueID,  runID, sampleID,coverageID, vcallerID, assayID, instrumentID, environmentID, status from $table where queueID='$queueID';"

while read -r queueID  runID  sampleID coverageID vcallerID assayID instrumentID environmentID status;
do
	echo "processing job-----> $queueID"

	# echo " job is $queueID , $runID, $sampleID, $coverageID, $vcallerID,  $assayID, $instrumentID, $environmentID, $status"


	# first update sampleAnalysisQueue table and set status of this queue to 1 i.e started processing
	updatestatement="UPDATE $table SET status=1 WHERE queueID = $queueID;"
	mysql --user="$user" --password="$password" --database="$database" --execute="$updatestatement"


	# second insert into pipelineStatus table to update status as started
	insertstatement="INSERT INTO pipelineStatus (queueID, plStatus, timeUpdated) VALUES ('$queueID','started',now());"
	mysql --user="$user" --password="$password" --database="$database" --execute="$insertstatement"

	# run pipeline

	if [ "$instrumentID" == "proton" ] || [ "$instrumentID" == "pgm" ]
	then
  	echo "proton or pgm"
  	bash /home/pipelines/master/shell/ionPipelineInterface.sh -r $runID -s $sampleID -c $coverageID -v $vcallerID -a $assayID -i $instrumentID -e $environmentID -q $queueID
	elif [ "$instrumentID" == "nextseq" ] || [ "$instrumentID" == "miseq" ]
	then
  	echo "nextseq or miseq"
  	bash /home/pipelines/master/shell/illuminaPipelineInterface.sh -r $runID -s $sampleID -a $assayID -i $instrumentID -e $environmentID -q $queueID
	fi

done < <(mysql --user="$user" --password="$password" --database="$database" --execute="$job_statement" -N)

DATE=`date '+%Y-%m-%d %H:%M:%S'`
echo "pipeline thread completed for job $queueID - $DATE" >> /home/cron_log.txt

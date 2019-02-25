#!/usr/bin/env bash
SHELL=/bin/bash

if [ $# -eq 0 ]
then
	echo "Usage: pipelineThread.sh"
	echo "-q queueID"
	echo "-e environment"
	echo "-u user"
	echo "-p password"
	exit
fi




if test $# -gt 0
	then
	while getopts :q:e:u:p: opt
	do
	case $opt in
  q)
		queueID=$OPTARG
		;;
	e)
		environment=$OPTARG
		;;
	u)
	  user=$OPTARG
	  ;;
	p)
	  password=$OPTARG
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

echo "Time: $(date -Iseconds) - INFO :  -q "$queueID" -e "$environment" -u "$user" -p "$password" "


echo "Time: $(date -Iseconds) - INFO : running pipelineThread.sh"
################################################################################
# running pipeline from sampleAnalysisQueue
################################################################################

job_statement="select pipelineQueue.queueID, samples.runID, samples.sampleName, samples.coverageID, samples.callerID, assays.assayName, instruments.instrumentName, pipelineQueue.status from pipelineQueue join samples on samples.sampleID=pipelineQueue.sampleID join assays on assays.assayID = samples.assayID join instruments on instruments.instrumentID = samples.instrumentID where pipelineQueue.queueID='$queueID';"

while read -r queueID  runID  sampleName coverageID callerID assay instrument status;
do
	echo "Time: $(date -Iseconds) - INFO : processing job -----> $queueID"
	echo "Time: $(date -Iseconds) - INFO : job is $queueID , $runID, $sampleName, $coverageID, $callerID,  $assay, $instrument, $environment, $status"

	# insert into pipelineStatus table to update status as started
	insertstatement="INSERT INTO pipelineStatus (queueID, plStatus, timeUpdated) VALUES ('$queueID','started',now());"
	mysql --user="$user" --password="$password" --database="$environment" --execute="$insertstatement"

	# run pipeline

	if [ "$instrument" == "proton" ] || [ "$instrument" == "pgm" ]
	then
		echo "Time: $(date -Iseconds) - INFO : running ionPipelineInterface.sh"
  	bash /var/pipelines_"$environment"/shell/ionPipelineInterface.sh -r $runID -s $sampleName -c $coverageID -v $callerID -a $assay -i $instrument -e $environment -q $queueID -u $user -p $password
	elif [ "$instrument" == "nextseq" ] || [ "$instrument" == "miseq" ]
	then
		echo "Time: $(date -Iseconds) - INFO : running illuminaPipelineInterface.sh"
  	bash /var/pipelines_"$environment"/shell/illuminaPipelineInterface.sh -r $runID -s $sampleName -a $assay -i $instrument -e $environment -q $queueID -u $user -p $password
	fi

done < <(sudo mysql --user="$user" --password="$password" --database="$environment" --execute="$job_statement" -N )

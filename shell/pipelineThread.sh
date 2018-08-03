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

echo " -q "$queueID" -e "$environment" -u "$user" -p "$password" "

################################################################################
# running pipeline from sampleAnalysisQueue
################################################################################

job_statement="select queueID,  runID, sampleID,coverageID, vcallerID, assayID, instrumentID, environmentID, status from sampleAnalysisQueue where queueID='$queueID';"

while read -r queueID  runID  sampleID coverageID vcallerID assayID instrumentID environmentID status;
do
	echo "processing job-----> $queueID"


	echo " job is $queueID , $runID, $sampleID, $coverageID, $vcallerID,  $assayID, $instrumentID, $environmentID, $status"
	echo "instrumentID is "$instrumentID" "

	# insert into pipelineStatus table to update status as started
	insertstatement="INSERT INTO pipelineStatus (queueID, plStatus, timeUpdated) VALUES ('$queueID','started',now());"
	mysql --user="$user" --password="$password" --database="$environment" --execute="$insertstatement"

	# run pipeline

	if [ "$instrumentID" == "proton" ] || [ "$instrumentID" == "pgm" ]
	then
  	echo "running -- instrument $instrumentID -- assay $assayID -- run $runID -- sample id is $sampleID"
		echo "running ionPipelineInterface.sh"
  	bash /home/pipelines/master/shell/ionPipelineInterface.sh -r $runID -s $sampleID -c $coverageID -v $vcallerID -a $assayID -i $instrumentID -e $environmentID -q $queueID -u $user -p $password
	elif [ "$instrumentID" == "nextseq" ] || [ "$instrumentID" == "miseq" ]
	then
  	echo "running -- instrument $instrumentID -- assay $assayID -- run $runID -- sample id is $sampleID"
		echo "instrumentID is "$instrumentID" "
		echo "running illuminaPipelineInterface.sh"
  	bash /home/pipelines/master/shell/illuminaPipelineInterface.sh -r $runID -s $sampleID -a $assayID -i $instrumentID -e $environmentID -q $queueID -u $user -p $password
	fi

done < <(sudo mysql --user="$user" --password="$password" --database="$environment" --execute="$job_statement" -N )

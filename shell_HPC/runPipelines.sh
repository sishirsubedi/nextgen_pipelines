#!/usr/bin/env bash
export SHELL=/usr/bin/bash

if [ $# -eq 0 ]
then
	echo "Usage: runPipelines.sh"
	echo "-e ENVIRONMENT"
  echo "-u USER"
  echo "-p PASSWORD"
	exit
fi

if test $# -gt 0
	then
	while getopts :e:u:p: opt
	do
	case $opt in
  e)
	ENVIRONMENT=$OPTARG
	;;
  u)
  USER=$OPTARG
  ;;
  p)
  PASSWORD=$OPTARG
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

HOME="/home/pipelines/pipelines_${ENVIRONMENT}/"
MYSQL_HOST="hhplabngsp01"

#remove previous empty log files
rm -f $(find ${HOME}cron_logs/ -name "*.log" -type f -size -375c)


date_=`date '+%Y-%m-%d %H:%M:%S'`
currentdate=$(echo $date_ | sed -e 's/ /_/g' -e 's/:/_/g' -e 's/-/_/g' )
touch  ${HOME}run_files/"$currentdate".txt
chmod 777  ${HOME}run_files/"$currentdate".txt

echo "Time: $(date -Iseconds) - INFO : Running pipelines"


################################################################################
# invoke nextseq
################################################################################

# nohup bash ${HOME}shell/runPipelines_nextseq.sh -e $ENVIRONMENT -u $USER -p $PASSWORD & > ${HOME}cron_logs/"nextseq_$currentdate".txt

################################################################################
# running non nextseq pipeline from sampleAnalysisQueue
################################################################################
defaultStatus=0
instrument='nextseq'
statement="select pipelineQueue.queueID from pipelineQueue join samples on samples.sampleID=pipelineQueue.sampleID \
					join assays on assays.assayID = samples.assayID join instruments on instruments.instrumentID = samples.instrumentID \
					where pipelineQueue.status='$defaultStatus'and instruments.instrumentName !='$instrument' order by pipelineQueue.queueID;"

echo "Time: $(date -Iseconds)- INFO - Queued MISEQ/PROTON jobs"
echo "--host="$MYSQL_HOST" --user="$USER" --password="$PASSWORD" --database="$ENVIRONMENT""
while  read -r queueID ;
do

   echo "Time: $(date -Iseconds)- INFO - MISEQ/PROTON queueID is $queueID"

   #update sampleAnalysisQueue table and set status of this queue to 1 i.e started processing
   updatestatement="UPDATE pipelineQueue SET status=1 WHERE queueID = $queueID;"
	 mysql --host="$MYSQL_HOST" --user="$USER" --password="$PASSWORD" --database="$ENVIRONMENT" --execute="$updatestatement"


   echo "$queueID" # >> ${HOME}run_files/"$currentdate".txt

done < <(mysql --host="$MYSQL_HOST" --user="$USER" --password="$PASSWORD" --database="$ENVIRONMENT" --execute="$statement" -N)


## run jobs in parallel
# /opt/parallel/bin/parallel --jobs ${HOME}run_files/jobfile \
#          -a ${HOME}run_files/"$currentdate".txt  "${HOME}shell/pipelineThread.sh -e $ENVIRONMENT -u $USER -p $PASSWORD -q "
#

# sudo rm -f ${HOME}run_files/"$currentdate".txt

echo "Time: $(date -Iseconds)- INFO - cron job finished running"

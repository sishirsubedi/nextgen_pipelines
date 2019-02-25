#!/usr/bin/env bash
export SHELL=/usr/bin/bash

if [ $# -eq 0 ]
then
	echo "Usage: runPipelines.sh"
	echo "-e environment"
  echo "-u user"
  echo "-p password"

	exit
fi

if test $# -gt 0
	then
	while getopts :e:u:p: opt
	do
	case $opt in
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

#remove previous empty log files
rm -f $(find /var/pipelines_"$environment"/cron_logs/ -name "*.log" -type f -size -375c)


date_=`date '+%Y-%m-%d %H:%M:%S'`
currentdate=$(echo $date_ | sed -e 's/ /_/g' -e 's/:/_/g' -e 's/-/_/g' )
echo "Time: $(date -Iseconds) - INFO : Running pipelines"
touch  /var/pipelines_"$environment"/run_files/"$currentdate".txt
sudo chmod 777  /var/pipelines_"$environment"/run_files/"$currentdate".txt


################################################################################
# invoke nextseq
################################################################################

nohup bash /var/pipelines_"$environment"/shell/runPipelines_nextseq.sh -e $environment -u $user -p $password & > /var/pipelines_"$environment"/cron_logs/"nextseq_$currentdate".txt

################################################################################
# running non nextseq pipeline from sampleAnalysisQueue
################################################################################
defaultStatus=0
instrument='nextseq'
statement="select pipelineQueue.queueID from pipelineQueue join samples on samples.sampleID=pipelineQueue.sampleID \
					join assays on assays.assayID = samples.assayID join instruments on instruments.instrumentID = samples.instrumentID \
					where pipelineQueue.status='$defaultStatus'and instruments.instrumentName !='$instrument' order by pipelineQueue.queueID;"

echo "Time: $(date -Iseconds)- INFO - Queued MISEQ/PROTON jobs"
while  read -r queueID ;
do

   echo "Time: $(date -Iseconds)- INFO - MISEQ/PROTON queueID is $queueID"

   #update sampleAnalysisQueue table and set status of this queue to 1 i.e started processing
   updatestatement="UPDATE pipelineQueue SET status=1 WHERE queueID = $queueID;"

   mysql --user="$user" --password="$password" --database="$environment" --execute="$updatestatement"
   echo "$queueID"  >> /var/pipelines_"$environment"/run_files/"$currentdate".txt

done < <(mysql --user="$user" --password="$password" --database="$environment" --execute="$statement" -N)


## run jobs in parallel
/opt/parallel/bin/parallel --jobs /var/pipelines_"$environment"/run_files/jobfile \
         -a /var/pipelines_"$environment"/run_files/"$currentdate".txt  "/var/pipelines_"$environment"/shell/pipelineThread.sh -e $environment -u $user -p $password -q "


sudo rm -f /var/pipelines_"$environment"/run_files/"$currentdate".txt

echo "Time: $(date -Iseconds)- INFO - cron job finished running"

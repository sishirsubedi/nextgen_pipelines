#!/usr/bin/env bash
export SHELL=/usr/bin/bash

if [ $# -eq 0 ]
then
	echo "Usage: pipelineThread.sh"
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


function updateStatus() {
user=$4
password=$5
database=$3
insertstatement="INSERT INTO pipelineStatus (queueID, plStatus, timeUpdated) VALUES ('$1','$2',now());"
mysql --user="$user" --password="$password" --database="$database" --execute="$insertstatement"
}


date_=`date '+%Y-%m-%d %H:%M:%S'`
currentdate=$(echo $date_ | sed -e 's/ /_/g' -e 's/:/_/g' -e 's/-/_/g' )
echo "Running nextseq pipelines cron job at - $currentdate "
touch  /var/pipelines_"$environment"/run_files/"nextseq_$currentdate".txt
sudo chmod 777  /var/pipelines_"$environment"/run_files/"nextseq_$currentdate".txt


################################################################################
# running nextseq pipeline from sampleAnalysisQueue
################################################################################
defaultStatus=0
instrument='nextseq'
nextseq_statement="select pipelineQueue.queueID from pipelineQueue join samples on samples.sampleID=pipelineQueue.sampleID join assays on assays.assayID = samples.assayID join instruments on instruments.instrumentID = samples.instrumentID where pipelineQueue.status='$defaultStatus'and instruments.instrumentName ='$instrument' order by pipelineQueue.queueID;"

echo "   ----- Queued NEXTSEQ jobs in this cronjob are --------"

while  read -r queueID  runID  sampleID instrumentID ;
do

	echo "queueID is $queueID , user is $user,  database is $environment "

	#update sampleAnalysisQueue table and set status of this queue to 1 i.e started processing
	updatestatement="UPDATE sampleAnalysisQueue SET status=1 WHERE queueID = $queueID;"
	mysql --user="$user" --password="$password" --database="$environment" --execute="$updatestatement"


	## convert bcl2fastq
	lockdir=/var/pipelines_"$environment"/run_files/nextseq.lock
	if mkdir "$lockdir" ## mkdir is atomic, file or file variable is not
	then

		echo "successfully acquired lock: $lockdir"
		echo "bcl2fastq running for  run - $runID ... "
		runFolder=$(ls -d /home/$instrumentID/*_"$runID"_*)
		updateStatus "$queueID" "bcl2fastq" "$environment" "$user"  "$password"
		bcl2fastq --no-lane-splitting --runfolder-dir $runFolder --output-dir $runFolder/out1  &&  echo " bcl2fastq completed for  run - $runID"
		sudo rm -rf "$lockdir"

	else

		echo "cannot acquire lock, giving up on $lockdir"
		echo "nextseq cron job did not complete at - $currentdate"
		updateStatus "$queueID" "ERROR:Second-bcl2fastq" "$environment" "$user"  "$password"
		exit 1

	fi

	echo "$queueID"  >> /var/pipelines_"$environment"/run_files/"nextseq_$currentdate".txt

done < <(mysql --user="$user" --password="$password" --database="$environment" --execute="$nextseq_statement" -N)



# run nextseq jobs
/opt/parallel/bin/parallel --jobs /var/pipelines_"$environment"/run_files/nextseq_jobfile \
				 --eta \
         -a /var/pipelines_"$environment"/run_files/"nextseq_$currentdate".txt  "/var/pipelines_"$environment"/shell/pipelineThread.sh -e $environment -u $user -p $password -q "


sudo rm -f /var/pipelines_"$environment"/run_files/"nextseq_$currentdate".txt

echo "nextseq cron job finished running at - $currentdate"

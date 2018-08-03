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
touch  /home/pipelines/master/run_files/"nextseq_$currentdate".txt
sudo chmod 777  /home/pipelines/master/run_files/"nextseq_$currentdate".txt


################################################################################
# running nextseq pipeline from sampleAnalysisQueue
################################################################################
defaultStatus=0
instrument='nextseq'
nextseq_statement="select queueID,runID,  sampleID, instrumentID from sampleAnalysisQueue where status='$defaultStatus'and instrumentID='$instrument' order by queueID;"

echo "   ----- Queued NEXTSEQ jobs in this cronjob are --------"

while  read -r queueID  runID  sampleID instrumentID ;
do

	echo "queueID is $queueID , user is $user,  database is $environment "

	#update sampleAnalysisQueue table and set status of this queue to 1 i.e started processing
	updatestatement="UPDATE sampleAnalysisQueue SET status=1 WHERE queueID = $queueID;"
	mysql --user="$user" --password="$password" --database="$environment" --execute="$updatestatement"

        ## convert bcl2fastq
	if [ ! -f /home/$instrumentID/*_"$runID"_*/out1/"$sampleID"*_R1_001.fastq.gz ]
	then

		echo " bcl2fastq running for  run - $runID ... "

		runFolder=$(ls -d /home/$instrumentID/*_"$runID"_*)

		updateStatus "$queueID" "bcl2fastq" "$environment" "$user"  "$password"

		bcl2fastq --no-lane-splitting --runfolder-dir $runFolder --output-dir $runFolder/out1  &&  echo " bcl2fastq completed for  run - $runID"

	fi

	echo "$queueID"  >> /home/pipelines/master/run_files/"nextseq_$currentdate".txt

done < <(mysql --user="$user" --password="$password" --database="$environment" --execute="$nextseq_statement" -N)



## run nextseq jobs

sudo parallel --jobs /home/pipelines/master/run_files/nextseq_jobfile \
         --load /home/pipelines/master/run_files/nextseq_loadfile \
				 --noswap \
				 --eta \
				 --memfree /home/pipelines/master/run_files/nextseq_memfile \
         -a /home/pipelines/master/run_files/"nextseq_$currentdate".txt  "/home/pipelines/master/shell/pipelineThread.sh -e $environment -u $user -p $password -q "


sudo rm -f /home/pipelines/master/run_files/"nextseq_$currentdate".txt

echo "nextseq cron job finished running at - $currentdate"

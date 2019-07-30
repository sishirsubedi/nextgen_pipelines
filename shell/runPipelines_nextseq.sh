#!/usr/bin/env bash
export SHELL=/usr/bin/bash

if [ $# -eq 0 ]
then
	echo "Usage: runPipelines_nextseq.sh"
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
rm -f $(find /var/pipelines_"$environment"/cron_logs/ -name "nextseq_*.txt" -type f -size -375c)


function updateStatus() {
user=$4
password=$5
database=$3
insertstatement="INSERT INTO pipelineStatus (queueID, plStatus, timeUpdated) VALUES ('$1','$2',now());"
mysql --user="$user" --password="$password" --database="$database" --execute="$insertstatement"
}


date_=`date '+%Y-%m-%d %H:%M:%S'`
currentdate=$(echo $date_ | sed -e 's/ /_/g' -e 's/:/_/g' -e 's/-/_/g' )
echo "Time: $(date -Iseconds) - INFO : Running nextseq pipelines "
touch  /var/pipelines_"$environment"/run_files/"nextseq_$currentdate".txt
sudo chmod 777  /var/pipelines_"$environment"/run_files/"nextseq_$currentdate".txt


################################################################################
# running nextseq pipeline from pipelineQueue
################################################################################
defaultStatus=0
instrument='nextseq'
nextseq_statement="select pipelineQueue.queueID , samples.runID from pipelineQueue join samples on samples.sampleID=pipelineQueue.sampleID join assays on assays.assayID = samples.assayID join instruments on instruments.instrumentID = samples.instrumentID where pipelineQueue.status='$defaultStatus'and instruments.instrumentName ='$instrument' order by pipelineQueue.queueID;"

echo "Time: $(date -Iseconds) - INFO : Queued NEXTSEQ jobs are"

ids=()


while  read -r queueID runID ;
do

	echo "Time: $(date -Iseconds) - INFO : Updating queue status for NEXTSEQ queueID - $queueID"

	#update sampleAnalysisQueue table and set status of this queue to 1 i.e started processing
	updatestatement="UPDATE pipelineQueue SET status=1 WHERE queueID = $queueID;"
	mysql --user="$user" --password="$password" --database="$environment" --execute="$updatestatement"

	echo "$queueID"  >> /var/pipelines_"$environment"/run_files/"nextseq_$currentdate".txt

	ids+=("$runID/$queueID")

done < <(mysql --user="$user" --password="$password" --database="$environment" --execute="$nextseq_statement" -N)


#uniqueRunIDS=($(printf "%s\n" "${runIDs[@]}" | sort -u))

################################################################################
# running bcl2fastq for unique runIDs
################################################################################

for idpair in "${ids[@]}"
do

	runID=$(echo $idpair | cut -d "/" -f 1)

	queueID=$(echo $idpair | cut -d "/" -f 2)

  fastqStatus=$(mysql --user="$user" --password="$password" --database="$environment" -se "select status from pipelineStatusBcl2Fastq where runID='$runID'")

  echo "Time: $(date -Iseconds) - INFO : checking bcl2fastq for runID - $runID and queueid $queueID and runid fastqStatus '$fastqStatus' --"

	if [ $fastqStatus == "0" ]
	then
    echo "Time: $(date -Iseconds) - INFO : running bcl2fastq for runID - $runID and queueid $queueID and runid fastqStatus '$fastqStatus' --"
		## convert bcl2fastq
		lockdir=/var/pipelines_"$environment"/run_files/nextseq.lock
		if mkdir "$lockdir" ## mkdir is atomic, file or file variable is not
		then

			updateStatus "$queueID" "bcl2fastq_running_now" "$environment" "$user"  "$password"

			echo "Time: $(date -Iseconds) - INFO : NEXTSEQ successfully acquired lock: $lockdir"
			echo "Time: $(date -Iseconds) - INFO : NEXTSEQ bcl2fastq running for  run - $runID queueID - $queueID ... "

			runFolder=$(ls -d /home/nextseq/*_"$runID"_*)

			sudo /usr/local/bin/bcl2fastq --no-lane-splitting --runfolder-dir $runFolder --output-dir $runFolder/out1 -r 5 -d 5 -p 5 -w 5 &

			wait $!

			updatestatement="UPDATE pipelineStatusBcl2Fastq SET status=1 WHERE runID = $runID;"
			mysql --user="$user" --password="$password" --database="$environment" --execute="$updatestatement"

			updateStatus "$queueID" "bcl2fastq_completed_now" "$environment" "$user"  "$password"

			sudo rm -rf "$lockdir"

		else

			echo "Time: $(date -Iseconds) - INFO : NEXTSEQ cannot acquire lock for bcl2fastq conversion $runID, giving up on $lockdir"
			echo "Time: $(date -Iseconds) - INFO : NEXTSEQ nextseq cron job DID NOT complete at - $currentdate"

			#update pipelineQueue table and set status of this queue to 0 i.e try again next time
			for idpair_re in "${ids[@]}"
			do

				runID_re=$(echo $idpair_re | cut -d "/" -f 1)
				queueID_re=$(echo $idpair_re | cut -d "/" -f 2)

				if [ $runID_re == $runID ]
				then

					updatestatement="UPDATE pipelineQueue SET pipelineQueue.status=0 WHERE queueID = $queueID_re;"
					mysql --user="$user" --password="$password" --database="$environment" --execute="$updatestatement"

					updateStatus "$queueID_re" "bcl2fastq_wait" "$environment" "$user"  "$password"
				fi

			done

			exit 1

		fi

	else

		echo "Time: $(date -Iseconds) - INFO : bcl2fastq_completed_past for runID - $runID and queueid $queueID and runid fastqStatus '$fastqStatus' --"
		updateStatus "$queueID" "bcl2fastq_completed_past" "$environment" "$user"  "$password"
		# this message means sample is already queued and waiting to run

fi
done


# run nextseq jobs
/opt/parallel/bin/parallel --jobs /var/pipelines_"$environment"/run_files/nextseq_jobfile \
         -a /var/pipelines_"$environment"/run_files/"nextseq_$currentdate".txt  "/var/pipelines_"$environment"/shell/pipelineThread.sh -e $environment -u $user -p $password -q "


sudo rm -f /var/pipelines_"$environment"/run_files/"nextseq_$currentdate".txt

echo "Time: $(date -Iseconds) - INFO : NEXTSEQ cron job finished running"

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

#remove previous empty log files
rm -f $(find /home/cron_logs/pipelines/ -name "*.log" -type f -size -250c)

date_=`date '+%Y-%m-%d %H:%M:%S'`
currentdate=$(echo $date_ | sed -e 's/ /_/g' -e 's/:/_/g' -e 's/-/_/g' )
echo "Running pipelines cron job at - $currentdate "
touch  /home/pipelines/master/run_files/"$currentdate".txt
sudo chmod 777  /home/pipelines/master/run_files/"$currentdate".txt


################################################################################
# invoke nextseq
################################################################################

nohup bash /home/pipelines/master/shell/runPipelines_nextseq.sh -e $environment -u $user -p $password &

sleep 2s
################################################################################
# running non nextseq pipeline from sampleAnalysisQueue
################################################################################
defaultStatus=0
instrument='nextseq'
statement="select queueID from sampleAnalysisQueue where status='$defaultStatus'and instrumentID !='$instrument' order by queueID;"
echo "   ----- Queued MISEQ/PROTON jobs in this cronjob are --------"
while  read -r queueID ;
do
   echo "queueID is $queueID , user is $user,  database is $environment "

   #update sampleAnalysisQueue table and set status of this queue to 1 i.e started processing
   updatestatement="UPDATE sampleAnalysisQueue SET status=1 WHERE queueID = $queueID;"
   mysql --user="$user" --password="$password" --database="$environment" --execute="$updatestatement"


   echo "$queueID"  >> /home/pipelines/master/run_files/"$currentdate".txt


done < <(mysql --user="$user" --password="$password" --database="$environment" --execute="$statement" -N)


sudo parallel --jobs /home/pipelines/master/run_files/jobfile \
         --load /home/pipelines/master/run_files/loadfile \
	       --noswap \
	       --eta \
	       --memfree /home/pipelines/master/run_files/memfile \
         -a /home/pipelines/master/run_files/"$currentdate".txt  "/home/pipelines/master/shell/pipelineThread.sh -e $environment -u $user -p $password -q "


sudo rm -f /home/pipelines/master/run_files/"$currentdate".txt


echo "cron job finished running at - $currentdate"

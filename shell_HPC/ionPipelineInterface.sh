#!/bin/bash

if [ $# -eq 0 ]
then
	echo "Usage: ionPipeline_Interface.sh"
	echo "-r runID: three digits, last number of the run"
	echo "-s sampleName"
	echo "-c coverageID: coverage analysis"
	echo "-v callerID: variant caller"
	echo "-a assay"
	echo "-i instrument"
	echo "-e environment"
	echo "-q queueID"
	echo "-u user"
	echo "-p password"
	exit
fi


if test $# -gt 0
	then
	while getopts :r:s:c:v:a:i:e:q:u:p: opt
	do
	case $opt in
	r)
		runID=$OPTARG
		;;
	s)
	  sampleName=$OPTARG
		;;
	c)
	  coverageID=$OPTARG
	  ;;
	v)
		callerID=$OPTARG
		;;
	a)
		assay=$OPTARG
		;;
	i)
    instrument=$OPTARG
	  ;;
  e)
    environment=$OPTARG
	  ;;
  q)
		queueID=$OPTARG
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

function log() {
 MESSAGE=$1
 TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
 SCRIPT="ionPipeline_Interface"
 echo " [ $TIMESTAMP ] [ $SCRIPT ] : $MESSAGE "
}

function updateStatus() {
user=$4
password=$5
database="ngs_$3"
insertstatement="INSERT INTO pipelineStatus (queueID, plStatus, timeUpdated) VALUES ('$1','$2',now());"
mysql --host="hhplabngsp01" --user="$user" --password="$password" --database="$database" --execute="$insertstatement"
}


variantFolder=$(ls -d /home/$instrument/*$runID/plugin_out/"$callerID")
ampliconFolder=$(ls -d /home/$instrument/*"$runID"/plugin_out/"$coverageID")
runFolder=$(ls -d /home/$instrument/*$runID)
runName=${runFolder##*/}
ENV_HOME="/home/environments/ngs_${environment}/"
SCRIPT_HOME="/home/pipelines/ngs_${environment}/shell/"

log " Running ionPipeline Interface for :
assay : $assay
instrument : $instrument
runID : $runID
sampleName : $sampleName
coverageID : $coverageID
callerID : $callerID
environment : $environment
queueID : $queueID "


if [ -z $runID ] || [ -z $sampleName ] || [ -z $coverageID ] || [ -z $callerID ] || [ -z $assay ] || [ -z $instrument ] || [ -z $environment ]
then
	echo "Error: Please input required parameters-"
	echo "-r runID: three digits, last number of the run"
	echo "-s sampleName"
	echo "-c coverageID: coverage analysis"
	echo "-v callerID: variant caller"
	echo "-a assay"
	echo "-i instrument"
	echo "-e environment"
	exit
fi

if [ ! -d /home/$instrument/*"$runID"/plugin_out/"$coverageID"/ ]
then
	echo "coverage folder for run ID:$runID; coverage ID:$coverageID not found"
    updateStatus "$queueID" "ERROR:Instrument_coverageID" "$environment" "$user"  "$password"
	exit
fi

if [ ! -d /home/$instrument/*"$runID"/plugin_out/"$callerID"/ ]
then
	echo "variant folder for run ID:$runID; variant caller ID:$callerID not found"
	  updateStatus "$queueID" "ERROR:Instrument_callerID" $environment "$user"  "$password"
	exit
fi


if [ ! -d ${ENV_HOME}${instrument}Analysis/$runName ]
then
	mkdir ${ENV_HOME}${instrument}Analysis/$runName
fi
chmod 775 ${ENV_HOME}${instrument}Analysis/$runName

if [ ! -d ${ENV_HOME}${instrument}Analysis/$runName/$sampleName ]
then
	mkdir ${ENV_HOME}${instrument}Analysis/$runName/$sampleName
fi
chmod 777 ${ENV_HOME}${instrument}Analysis/$runName/$sampleName


exec >  >(tee -a ${ENV_HOME}${instrument}Analysis/$runName/$sampleName/process.log)
exec 2> >(tee -a ${ENV_HOME}${instrument}Analysis/$runName/$sampleName/process.log >&2)

# running ionPipeline scripts based on assay and instrument
if [ $assay == "neuro" ]
then

	bash ${SCRIPT_HOME}ionPipeline.sh -r $runID -s $sampleName -c $coverageID -v $callerID -i $instrument  -e /home/doc/ref/neuralRef/excludedAmplicon.txt -a /home/doc/ref/neuralRef/IAD87786_179_Designed.excluded.bed -n $environment -q $queueID -u $user -p $password
	exit

elif [ $assay == "gene50" ]
then

	bash ${SCRIPT_HOME}ionPipeline.sh -r $runID -s $sampleName -c $coverageID -v $callerID -i $instrument -n $environment -q $queueID -u $user -p $password
	exit

else

	echo "Error: Failed ionPipeline for:
	assay : $assay
	instrument : $instrument
	runID : $runID
	sampleName : $sampleName
	coverageID : $coverageID
	callerID : $callerID
	environment : $environment
	queueID : $queueID "
	echo "Not valid assay - $assay and instrument - $instrument. Process Terminated."
  updateStatus "$queueID" "ERROR:assay_inst" $environment "$user"  "$password"
	exit
fi

done

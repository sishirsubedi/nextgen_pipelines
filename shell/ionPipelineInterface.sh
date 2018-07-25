#!/bin/bash

if [ $# -eq 0 ]
then
	echo "Usage: ionPipeline_Interface.sh"
	echo "-r runID: three digits, last number of the run"
	echo "-s sampleID"
	echo "-c coverageID: coverage analysis"
	echo "-v callerID: variant caller"
	echo "-a assayID"
	echo "-i instrumentID"
	echo "-e environmentID"
	echo "-q queueID"
	exit
fi


if test $# -gt 0
	then
	while getopts :r:s:c:v:a:i:e:q: opt
	do
	case $opt in
	r)
		runID=$OPTARG
		;;
	s)
	  sampleID=$OPTARG
		;;
	c)
	  coverageID=$OPTARG
	  ;;
	v)
		callerID=$OPTARG
		;;
	a)
		assayID=$OPTARG
		;;
	i)
    instrumentID=$OPTARG
	  ;;
  e)
    environmentID=$OPTARG
	  ;;
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




function updateStatus() {

user=hhadmin
password=ngs3127
database=$3
insertstatement="INSERT INTO pipelineStatus (queueID, plStatus, timeUpdated) VALUES ('$1','$2',now());"
mysql --user="$user" --password="$password" --database="$database" --execute="$insertstatement"
}





if [ -z $runID ] || [ -z $sampleID ] || [ -z $coverageID ] || [ -z $callerID ] || [ -z $assayID ] || [ -z $instrumentID ] || [ -z $environmentID ]
then
	echo "Error: Please input required parameters-"
	echo "-r runID: three digits, last number of the run"
	echo "-s sampleID"
	echo "-c coverageID: coverage analysis"
	echo "-v callerID: variant caller"
	echo "-a assayID"
	echo "-i instrumentID"
	echo "-e environmentID"

   #ERROR:parameters

	exit
fi

if [ ! -d /home/$instrumentID/*"$runID"/plugin_out/coverageAnalysis_out."$coverageID"/ ]
then
	echo "coverage folder for run ID:$runID; coverage ID:$coverageID not found"
    updateStatus "$queueID" "ERROR:Instrument_coverageID" "$environmentID"
	exit
fi

if [ ! -d /home/$instrumentID/*"$runID"/plugin_out/variantCaller_out."$callerID"/ ]
then
	echo "variant folder for run ID:$runID; variant caller ID:$callerID not found"
	  updateStatus "$queueID" "ERROR:Instrument_callerID" $environmentID
	exit
fi


variantFolder=$(ls -d /home/$instrumentID/*$runID/plugin_out/variantCaller_out."$callerID")
ampliconFolder=$(ls -d /home/$instrumentID/*"$runID"/plugin_out/coverageAnalysis_out."$coverageID")
runFolder=$(ls -d /home/$instrumentID/*$runID)
runName=${runFolder##*/}

if [ ! -d /home/environments/$environmentID/"$instrumentID"Analysis/$runName ]
then
	mkdir /home/environments/$environmentID/"$instrumentID"Analysis/$runName
fi
chmod 777 /home/environments/$environmentID/"$instrumentID"Analysis/$runName

if [ ! -d /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID ]
then
	mkdir /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID
fi
chmod 777 /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID


exec >  >(tee -a /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/process.log)
exec 2> >(tee -a /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/process.log >&2)

# running ionPipeline scripts based on assay and instrument
if [ $assayID == "neuro" ] && [ $instrumentID == "proton" ]
then
	echo "Running ionPipeline for :
	assayID : $assayID
	instrumentID : $instrumentID
	runID : $runID
	sampleID : $sampleID
	coverageID : $coverageID
	callerID : $callerID
	environmentID : $environmentID
	queueID : $queueID "

	bash /home/pipelines/master/shell/ionPipeline.sh -r $runID -s $sampleID -c $coverageID -v $callerID -i $instrumentID  -e /home/doc/ref/neuralRef/excludedAmplicon.txt -a /home/doc/ref/neuralRef/IAD87786_179_Designed.excluded.bed -n $environmentID -q $queueID

	exit

elif [ $assayID == "gene50" ] && ( [ $instrumentID == "proton" ] || [ $instrumentID == "pgm" ] )
then
	echo "Running ionPipeline for:
	assayID : $assayID
	instrumentID : $instrumentID
	runID : $runID
	sampleID : $sampleID
	coverageID : $coverageID
	callerID : $callerID
	environmentID : $environmentID
  queueID : $queueID "

	bash /home/pipelines/master/shell/ionPipeline.sh -r $runID -s $sampleID -c $coverageID -v $callerID -i $instrumentID -n $environmentID -q $queueID
	exit

else

	echo "Error: Failed ionPipeline for:
	assayID : $assayID
	instrumentID : $instrumentID
	runID : $runID
	sampleID : $sampleID
	coverageID : $coverageID
	callerID : $callerID
	environmentID : $environmentID
	queueID : $queueID "
	echo "Not valid assay - $assayID and instrument - $instrumentID. Process Terminated."
    updateStatus "$queueID" "ERROR:assay_inst" $environmentID
	exit
fi

done

#!/bin/bash

if [ $# -eq 0 ]
then
	echo "Usage: illuminaPipeline_Interface.sh"
	echo "-r runID: three digits, last number of the run"
	echo "-s sampleID"
	echo "-a assayID"
	echo "-i instrumentID"
	echo "-e environmentID"
	echo "-q queueID"
	exit
fi

if test $# -gt 0
	then
	while getopts :r:s:a:i:e:q: opt
	do
	case $opt in
	r)
		runID=$OPTARG
		;;
	s)
	  sampleID=$OPTARG
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

if [ -z $runID ] || [ -z $sampleID ] || [ -z $assayID ] || [ -z $instrumentID ] || [ -z $environmentID ]
then
	echo "Error: Please input required parameters-"
	echo "-r runID: three digits, last number of the run"
	echo "-s sampleID"
	echo "-a assayID"
	echo "-i instrumentID"
	echo "-e environmentID"
	exit
fi

if [ $assayID == "heme" ] && [ $instrumentID == "miseq" ]
then

  runFolder=$(ls -d /home/$instrumentID/*_"$runID"_*/Data/Intensities/BaseCalls/Alignment/)
	post=${runFolder##/home/$instrumentID/}
	runName=${post%%/Data/Intensities/BaseCalls/Alignment*}

	echo "runfolder $runFolder"
	echo "post $post"
	echo "runName $runName"


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


	echo "Running illuminaPipeline for :
	assayID : $assayID
	instrumentID : $instrumentID
	runID : $runID
	sampleID : $sampleID
	environmentID : $environmentID"

	bash /home/pipelines/master/shell/illuminaPipeline.sh -r $runID -s $sampleID -i $instrumentID  -e /doc/ref/Heme/excludedAmplicons.txt -a /doc/ref/Heme/trusight-myeloid-amplicon-track.excluded.bed -n $environmentID -q $queueID
	exit

elif [ $assayID == "heme" ] && [ $instrumentID == "nextseq" ]
then

	runFolder=$(ls -d /home/$instrumentID/*_"$runID"_*)      #eg: /home/nextseq/150807_NS500761_0011_AH3TTJAFXX
	echo $runFolder
	runName=${runFolder##/home/$instrumentID/}				#eg: 150807_NS500761_0011_AH3TTJAFXX
	echo $runName

	if [ ! -d /home/environments/$environmentID/"$instrumentID"_heme/$runName ]
	then
		mkdir /home/environments/$environmentID/"$instrumentID"_heme/$runName
	fi
	chmod 777 /home/environments/$environmentID/"$instrumentID"_heme/$runName

	if [ ! -d /home/environments/$environmentID/"$instrumentID"_heme/$runName/$sampleID ]
	then
		mkdir /home/environments/$environmentID/"$instrumentID"_heme/$runName/$sampleID
	fi
	chmod 777 /home/environments/$environmentID/"$instrumentID"_heme/$runName/$sampleID

	exec >  >(tee -a /home/environments/$environmentID/"$instrumentID"_heme/$runName/$sampleID/process.log)
	exec 2> >(tee -a /home/environments/$environmentID/"$instrumentID"_heme/$runName/$sampleID/process.log >&2)

	echo "Running illuminaPipeline for :
	assayID : $assayID
	instrumentID : $instrumentID
	runID : $runID
	sampleID : $sampleID
	environmentID : $environmentID"


  bash /home/pipelines/master/shell/updatepipelineStatus.sh -q $queueID -s bcl2fastq

	if [ ! -f /home/$instrumentID/*_"$runID"_*/out1/"$sampleID"*_R1_001.fastq.gz ]
	then
		echo "Fastq files not found. Converting bcl to fastq"
		bcl2fastq --no-lane-splitting --runfolder-dir $runFolder --output-dir $runFolder/out1
	fi

  bash /home/pipelines/master/shell/updatepipelineStatus.sh -q $queueID -s varscanPE

	##Aligning fastq files
	echo "Aligning fastq files"
	file=$runFolder/out1/"$sampleID"*_R1_001.fastq.gz
	fastq1=$file
	fastq2=${file/_R1_/_R2_}      #repleace "R1" with "R2"
  echo $file
	echo $fastq2

	bash /home/pipelines/master/shell/VarScanPipelinePE.sh -p $fastq1 -q $fastq2 -o /home/environments/$environmentID/"$instrumentID"_heme/$runName/$sampleID/


	# ##Variant calling
	echo "Starting variant calling"
	bash /home/pipelines/master/shell/illuminaPipeline.sh -r $runID -s $sampleID -i $instrumentID  -e /doc/ref/Heme/excludedAmplicons.txt -a /doc/ref/Heme/trusight-myeloid-amplicon-track.excluded.bed -n $environmentID -q $queueID

  exit

else
		echo "Error: Failed ionPipeline for:
		assayID : $assayID
		instrumentID : $instrumentID
		runID : $runID
		sampleID : $sampleID
		environmentID : $environmentID "
		echo "Not valid assay - $assayID and instrument - $instrumentID. Process Terminated."
		exit
fi

done

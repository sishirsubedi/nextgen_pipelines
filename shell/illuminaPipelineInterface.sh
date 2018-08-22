#!/bin/bash

if [ $# -eq 0 ]
then
	echo "Usage: illuminaPipeline_Interface.sh"
	echo "-r runID: three digits, last number of the run"
	echo "-s sampleName"
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
	while getopts :r:s:a:i:e:q:u:p: opt
	do
	case $opt in
	r)
		runID=$OPTARG
		;;
	s)
	  sampleName=$OPTARG
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


function updateStatus() {

user=$4
password=$5
database=$3
insertstatement="INSERT INTO pipelineStatus (queueID, plStatus, timeUpdated) VALUES ('$1','$2',now());"
mysql --user="$user" --password="$password" --database="$database" --execute="$insertstatement"
}


if [ -z $runID ] || [ -z $sampleName ] || [ -z $assay ] || [ -z $instrument ] || [ -z $environment ]
then
	echo "Error: Please input required parameters-"
	echo "-r runID: three digits, last number of the run"
	echo "-s sampleName"
	echo "-a assay"
	echo "-i instrument"
	echo "-e environment"
	  updateStatus "$queueID" "ERROR:parameters" "$environment"  "$user"  "$password"
	exit
fi

if [ $assay == "heme" ] && [ $instrument == "miseq" ]
then

  runFolder=$(ls -d /home/$instrument/*_"$runID"_*/Data/Intensities/BaseCalls/Alignment/)
	post=${runFolder##/home/$instrument/}
	runName=${post%%/Data/Intensities/BaseCalls/Alignment*}

	echo "runfolder $runFolder"
	echo "post $post"
	echo "runName $runName"


	if [ ! -d /home/environments/$environment/"$instrument"Analysis/$runName ]
	then
		mkdir /home/environments/$environment/"$instrument"Analysis/$runName
	fi
	chmod 777 /home/environments/$environment/"$instrument"Analysis/$runName

	if [ ! -d /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName ]
	then
		mkdir /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName
	fi
	chmod 777 /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName

	exec >  >(tee -a /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/process.log)
	exec 2> >(tee -a /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/process.log >&2)


	echo "Running illuminaPipeline for :
	assay : $assay
	instrument : $instrument
	runID : $runID
	sampleName : $sampleName
	environment : $environment"

	bash /var/pipelines_"$environment"/shell/illuminaPipeline.sh -r $runID -s $sampleName -i $instrument  -e /doc/ref/Heme/excludedAmplicons.txt -a /doc/ref/Heme/trusight-myeloid-amplicon-track.excluded.bed -n $environment -q $queueID -u $user -p $password
	exit

elif [ $assay == "heme" ] && [ $instrument == "nextseq" ]
then

	runFolder=$(ls -d /home/$instrument/*_"$runID"_*)      #eg: /home/nextseq/150807_NS500761_0011_AH3TTJAFXX
	echo $runFolder
	runName=${runFolder##/home/$instrument/}				#eg: 150807_NS500761_0011_AH3TTJAFXX
	echo $runName

	if [ ! -d /home/environments/$environment/"$instrument"Analysis/$runName ]
	then
		mkdir /home/environments/$environment/"$instrument"Analysis/$runName
	fi
	chmod 777 /home/environments/$environment/"$instrument"Analysis/$runName


	if [ ! -d /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName ]
	then
		mkdir /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName
	fi
	chmod 777 /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName


	if [ ! -d /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantCaller ]
	then
		mkdir /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantCaller
	fi
	chmod 777 /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantCaller

	exec >  >(tee -a /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantCaller/process.log)
	exec 2> >(tee -a /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantCaller/process.log >&2)

	echo "Running illuminaPipeline for :
	assay : $assay
	instrument : $instrument
	runID : $runID
	sampleName : $sampleName
	environment : $environment"


	if [ ! -f /home/$instrument/*_"$runID"_*/out1/"$sampleName"*_R1_001.fastq.gz ]
	then
		echo "Fastq files not found"
		updateStatus "$queueID" "ERROR:fastqNotFound" "$environment" "$user"  "$password"
		exit
	fi

  updateStatus "$queueID" "varscanPE" "$environment" "$user"  "$password"

	##Aligning fastq files
	echo "Aligning fastq files"
	file=$runFolder/out1/"$sampleName"*_R1_001.fastq.gz
	fastq1=$file
	fastq2=${file/_R1_/_R2_}      #repleace "R1" with "R2"
  echo $file
	echo $fastq2

	bash /var/pipelines_"$environment"/shell/VarScanPipelinePE.sh -p $fastq1 -q $fastq2 -o /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantCaller


	# ##Variant calling
	echo "Starting variant calling"
	bash /var/pipelines_"$environment"/shell/illuminaPipeline.sh -r $runID -s $sampleName -i $instrument  -e /doc/ref/Heme/excludedAmplicons.txt -a /doc/ref/Heme/trusight-myeloid-amplicon-track.excluded.bed -n $environment -q $queueID -u $user -p $password

  exit

elif [ $assay == "exome" ] && [ $instrument == "nextseq" ]
then

	runFolder=$(ls -d /home/$instrument/*_"$runID"_*)
	echo $runFolder
	runName=${runFolder##/home/$instrument/}
	echo $runName

	if [ ! -d /home/environments/$environment/"$instrument"Analysis/$runName ]
	then
		mkdir /home/environments/$environment/"$instrument"Analysis/$runName
	fi
	chmod 777 /home/environments/$environment/"$instrument"Analysis/$runName


	if [ ! -d /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName ]
	then
		mkdir /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName
	fi
	chmod 777 /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName


	if [ ! -d /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantCaller ]
	then
		mkdir /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantCaller
	fi
	chmod 777 /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantCaller

	exec >  >(tee -a /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantCaller/process.log)
	exec 2> >(tee -a /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantCaller/process.log >&2)

	echo "Running illuminaPipeline for :
	assay : $assay
	instrument : $instrument
	runID : $runID
	sampleName : $sampleName
	environment : $environment"


	if [ ! -f /home/$instrument/*_"$runID"_*/out1/"$sampleName"*_R1_001.fastq.gz ]
	then
		echo "Fastq files not found"
		updateStatus "$queueID" "ERROR:fastqNotFound" "$environment" "$user"  "$password"
		exit
	fi

	updateStatus "$queueID" "varscanPE" "$environment" "$user"  "$password"

	##Aligning fastq files
	echo "Aligning fastq files"
	file=$runFolder/out1/"$sampleName"*_R1_001.fastq.gz
	fastq1=$file
	fastq2=${file/_R1_/_R2_}      #repleace "R1" with "R2"
  echo $file
	echo $fastq2

	bash /var/pipelines_"$environment"/shell/gatkPipeline.sh -p $fastq1 -q $fastq2 -o /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantCaller


	exit

else
		echo "Error: Failed ionPipeline for:
		assay : $assay
		instrument : $instrument
		runID : $runID
		sampleName : $sampleName
		environment : $environment "
		echo "Not valid assay - $assay and instrument - $instrument. Process Terminated."
		exit
fi

done

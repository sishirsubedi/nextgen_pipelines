#!/bin/sh
##check number of parameters##
if [ $# -eq 0 ]
then
echo "Usage: blastJob.sh"
echo "-r runID: third field in Run folder name"
echo "-n blastNumber: number of reads to blast ie. 500000, 500000, Default:0(full blast)"
echo "-s SampleID: The sample name typed into your SampleSheet, seperated by , ; If omitted, process all samples in the run" 
echo "-f input folder: Where the fastq files are, incompatible with -r and -s"
echo "-d database to blast against: default: micro(bacteria and fungi) Options: nt, myco, HN"
exit
fi

##check if current dir contains files with output file names##
if [ -f fileList.txt ] || [ -f outDirList.txt ] || [ -f sampleList.txt ] || [ -f torque ]
then
echo "Error: Please make sure no files named fileList.txt outDirList.txt sampleList.txt torque in current folder"
exit
fi

if test $# -gt 0
	then
	while getopts :r:n:s:d:f: opt
	do
	case $opt in
	
	n)
		blastNumber=$OPTARG
		;;
	r)
		runID=$OPTARG
		;;
	s)
		sampleID=$OPTARG
		;;
	d)
		dataBase=$OPTARG
		;;
	f)
		folder=$OPTARG
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
	

	
##check if all required arguments are supplied##
if [ -z $folder ] && [ -z $runID ] && [ -z $sampleID ]
then
echo "Error: You need to specify either a folder where all fastq files are -f or both Run ID -r and Sample IDs -s"
exit
elif [ ! -z $folder ] && [ ! -z $runID ]
then
echo "Error: please specify either a folder -f  or Run ID -r and samples -s , not both"
exit
elif [ ! -z $folder ] && [ ! -z $sampleID ]
then
echo "Error: please specify either a folder -f or Run ID -r and samples -s, not both"
exit
elif [ -z $runID ] && [ ! -z $sampleID ]
then
echo "Error: please specify a run ID -r"
exit
elif [ ! -z $runID ] && [ -z $sampleID ]
then
echo "Will run blast on all samples for run $runID"
fi 

###assign default database and blastNumber##
dataBase=${dataBase:-micro}
blastNumber=${blastNumber:-0}

##not the correct folder name##
if [ ! -z $runID ] && [ ! -d /home/miseq/*_"$runID"_*/ ]
then
echo "Error: run folder starting with $dir not found; nothing was done"
exit
fi	


#find the dir where the fastqs are
if [ ! -z $runID ]
then
fromDir=$(ls -d /home/miseq/*_"$runID"_*/Data/Intensities/BaseCalls)
else
fromDir=$folder
fi

####get fastq file info###
if [ ! -z $sampleID ]
then
IFS=,
sampleArray=($sampleID)
else
n=0
for file in $fromDir/*.fastq.gz
do
fileName=${file##*/}
sample=${fileName%%_*}
sampleArray[$n]="$sample"
n=$((n+1))
done
fi

##construct output directories
currentDir=$(pwd)
n=0
for key in "${!sampleArray[@]}"
do
sample=${sampleArray[$key]} 
echo $sample >> sampleList.txt
mkdir $sample 
##make file list##
echo "$currentDir/$sample" >> outDirList.txt
ls $fromDir/"$sample"_*.fastq.gz >> fileList.txt
n=$((n+1))
done


if [[ $n -ge 8 ]]
then
nodes=4
n=8
else
nodes=$((n/2))
nodes=$((nodes+1))
fi
	
###construct the torque job script header####
echo "#!/bin/bash" >> torque
echo -e "\t#PBS -l walltime=99:00:00" >> torque
echo -e "\t#PBS -l nice=19" >> torque
echo -e "\t#PBS -q default" >> torque
echo -e "\t#PBS -l nodes=$nodes:ppn=12" >> torque


##write script###
echo "cd $currentDir" >> torque
echo 'cat $PBS_NODEFILE > nodes.txt' >> torque
echo "/opt/parallel/parallel -j 2 --sshloginfile nodes.txt -a fileList.txt -a outDirList.txt -a sampleList.txt --xapply 'sh /home/niyunyun/code/shell/blastSingle.sh -b $dataBase -f {1} -n $blastNumber -o {2}  > {2}/{3}.log 2>&1'" >> torque
echo "qsub torque"



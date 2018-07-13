if test $# -gt 0
	then
	while getopts :i:o:r:n:p: opt
	do
	case $opt in
	i)
		inDir=$OPTARG
		;;
	o)
		outDir=$OPTARG
		;;
	r)
		reference=$OPTARG
		;;
	n)
		nodes=$OPTARG
		;;
	p)
		ppn=$OPTARG
		;;
	:)
		echo "Option -$OPTARG requires an argument."
		;;
	\?)
		echo "Invalid option: -$OPTARG"
	esac
	done
	shift $((OPTIND-1))
	
if [ "$nodes" -gt 3 ] || [ "$ppn" -gt 9 ] 
then
	echo "maximum number nodes:3 and processes per node:9"
else
	if [ -e torque ] || [ -e fileList.txt ] || [ -e outDirList.txt ] || [ -e sampleList.txt ] || [ -e nodes.txt ]
	then
		echo "Aborted: please make sure that the working directory does not contain files with the following names:"
		echo "[torque] [fileList.txt] [outDirList.txt] [sampleList.txt] [nodes.txt]"
		exit
	fi
	totalProcesses=`expr $nodes \\* $ppn `
	###construct the torque job script header####
	echo "#!/bin/bash" >> torque
	echo -e "\t#PBS -l walltime=99:00:00" >> torque
	echo -e "\t#PBS -l nice=19" >> torque
	echo -e "\t#PBS -q default" >> torque
	echo -e "\t#PBS -l nodes=$nodes:ppn=$ppn" >> torque
	###construct the script###
	##make file list##
	ls $inDir/*.fastq.gz > fileList.txt
	##make sample list##
	for file in $inDir/*.fastq.gz
	do
	fileName=${file##*/}
	sample=${fileName%%.fastq.gz}
	echo $sample >> sampleList.txt
	mkdir $outDir/$sample
	echo $outDir/$sample >> outDirList.txt
	done
	
	##write script###
	echo "cd $outDir" >> torque
	echo 'cat $PBS_NODEFILE > nodes.txt' >> torque
	echo "/opt/parallel/parallel -j $totalProcesses --sshloginfile nodes.txt -a fileList.txt -a outDirList.txt -a sampleList.txt --xapply 'sh /home/niyunyun/code/shell/outBreak.sh {1} {2} $reference  > {2}/{3}.log 2>&1'" >> torque
	echo "qsub torque"
fi

else
	echo "Usage: sh VarScanJobPE.sh -i [input folder with fastq files] -o [output folder] -r [reference file] -n [number of nodes] (max:3) -p [number of processes per node] (max:9)"
fi	
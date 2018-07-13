#!/bin/sh
if [ $# -eq 0 ]
then
echo "Usage: readCountBatch.sh"
echo "-l hotspot file, first 4 columns: chr, pos, ref, alt"
echo "-d output directory"
echo "-i input directory with one or more bam files"
echo " output file name will be parsed from input bam file names"
exit
fi


if test $# -gt 0
	then
	while getopts :l:d:i: opt
	do
	case $opt in
	
	l)
		hotspot=$OPTARG
		;;
	d)
		outputDir=$OPTARG
		;;
	i)
		inputDir=$OPTARG
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


for file in $inputDir/*.bam
do
echo "processing $file"
fileName=${file##*/}
sampleName=${fileName%%.bam}
bash ~/code/shell/readCount.sh -l $hotspot -d $outputDir -o $sampleName -b $file
done

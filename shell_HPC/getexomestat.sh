#!/bin/bash

if [ $# -eq 0 ]
then
	echo "Please provide the following information to run get exome stat"
	echo "-d Run directory in nextseq folder"
	echo "-s sample"
	echo "Example:"
	echo "bash getexomestat.sh -d "190423_NS500761_0346_AHLVWTBGX9"  -s  "Exome14-N_S6"  "
	exit
fi

while getopts :d:s: opt; do
				case $opt in
        d)
  			DIR=$OPTARG
  			;;
			  s)
				SAMPLE=$OPTARG
				;;
			  :)
				echo "Option -$OPTARG requires an argument."
				;;
			  \?)
				echo "Invalid option: -$OPTARG"
		esac
done
	  shift $((OPTIND -1))

echo " exome job submitted for"
echo " folder - $DIR"
echo " sample - $SAMPLE"


/opt/python3/bin/python3 /home/hhadmin/scripts/bioinfoTools/02_bamQC/10_getSeqStat.py  "/home/environments/ngs_test/exomeAnalysis/${DIR}/Single/${SAMPLE}/Alignment/"   "${SAMPLE}"

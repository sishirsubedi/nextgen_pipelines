#!/bin/bash

if [ $# -eq 0 ]
then
	echo "Please provide the following information to run exome analysis"
	echo "-d Run directory in nextseq folder"
	echo "-n normal sample"
  echo "-t tumor sample"
	echo "-v Variant Calling"
	echo "-e environment"
	echo "-u USER"
	echo "-p PASSWORD"
	echo "RULE:
-v=NO   &&  -t=NO then alignment only
-v=YES  &&  -t=NO then unpaired variant caller
-v=YES  &&  -t=ExomeID then  paired variant caller"
	echo "Example:
bash /home/pipelines/ngs_test/shell/runexome.sh -d "190423_NS500761_0346_AHLVWTBGX9"  -n  "Exome14-N_S6"  -t "Exome14-T_S7"  -v YES -e test"
	exit
fi

while getopts :d:n:t:v:e:u:p: opt; do
				case $opt in
        d)
  			DIR=$OPTARG
  			;;
			  n)
				NORMAL=$OPTARG
				;;
			  t)
			  TUMOR=$OPTARG
			  ;;
				v)
			  VC=$OPTARG
			  ;;
			  e)
				ENVIRONMENT=$OPTARG
				;;
				u)
				USER=$OPTARG
				;;
				p)
				PASSWORD=$OPTARG
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
echo " normal - $NORMAL"
echo " tumor - $TUMOR"
echo " variantcalling - $VC"
echo " environment - $ENVIRONMENT"
echo " environment - $USER"

qsub -k eo -F "-d$DIR -n$NORMAL -t$TUMOR -v$VC -e$ENVIRONMENT -u$USER -p$PASSWORD" /home/pipelines/ngs_test/shell/exomeInterface.sh
#!/bin/bash
#PBS -l nodes=n001.cluster.com
#PBS -l walltime=99:00:00
#PBS -q default
#PBS -k eo


while getopts :c:e:u:p: opt; do
				case $opt in
        c)
  			FILE_ID=$OPTARG
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

INSTRUMENT="proton"
ASSAY="gene50"
HOME="/home/pipelines/ngs_${ENVIRONMENT}/"
HOME_RUNDIR="${HOME}run_files/"
HOME_SHELLDIR="${HOME}shell/"
DB="ngs_${ENVIRONMENT}"


function log() {
 MESSAGE=$1
 TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
 SCRIPT="submitJob_proton_gene50"
 echo " [ $TIMESTAMP ] [ $SCRIPT ] : $MESSAGE "
}


log "Running from QSUB:
ENVIRONMENT - $ENVIRONMENT
INSTRUMENT - $INSTRUMENT
ASSAY - Gene50
FILE_ID - $FILE_ID
SAMPLES -------------------------------------------------"
cat ${HOME_RUNDIR}proton_gene50_${FILE_ID}.samples| while IFS=';' read -ra line; do
echo "
RUN-ID : "${line[0]}"
SAMPLE-ID : "${line[1]}"
SAMPLE-NAME : "${line[2]}"
COVERAGE-ID : "${line[3]}"
CALLER-ID : "${line[4]}"
 "
done
echo " -------------------------------------------------"

cat ${HOME_RUNDIR}proton_gene50_${FILE_ID}.samples| while IFS=';' read -ra line; do
	QUEUEID="${line[0]}"
	RUNID="${line[1]}"
	SAMPLENAME="${line[2]}"
	COVERAGEID="${line[3]}"
	CALLERID="${line[4]}"
	bash ${HOME_SHELLDIR}ionPipelineInterface.sh -r $RUNID -s $SAMPLENAME -c $COVERAGEID -v $CALLERID -a $ASSAY -i $INSTRUMENT -e $ENVIRONMENT -q $QUEUEID -u $USER -p $PASSWORD
done

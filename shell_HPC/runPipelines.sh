##############################################################################
#
# Houston Methodist Hospital
# Molecular Diagnostic
#
#Script: NGS
#Description:
# This script is run by cronjob. It checks NGS database queue and submits jobs.
#
##############################################################################

#!/usr/bin/env bash

if [ $# -ne 6 ] ; then
	echo "Usage: runPipelines.sh";
	echo "Expected arguments:";
	echo "	-e ENVIRONMENT"
  echo "	-u USER"
  echo "	-p PASSWORD"
	exit
else
			while getopts :e:u:p: opt; do
				case $opt in
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
fi

################################################################################
# Initialize Variables
################################################################################

HOME="/home/pipelines/ngs_${ENVIRONMENT}/"
HOME_RUNDIR="${HOME}run_files/"
HOME_SHELLDIR="${HOME}shell/"
DB="ngs_${ENVIRONMENT}"

################################################################################
#
################################################################################

function log() {
 MESSAGE=$1
 TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
 SCRIPT="runPipelines"
 echo " [ $TIMESTAMP ] [ $SCRIPT ] : $MESSAGE "
}

function createSampleFiles() {
  current_dt=$1
	instrument=$2
	assay=$3
	touch  ${HOME_RUNDIR}${instrument}_${assay}_${current_dt}".samples"
	chmod 775  ${HOME_RUNDIR}${instrument}_${assay}_${current_dt}".samples"
}

function submitJob() {
  current_dt=$1
	instrument=$2
	assay=$3
  samfile="${HOME_RUNDIR}${instrument}_${assay}_${current_dt}".samples""
}


DEFAULT_STATUS=0
STATUS_QUERY_STATEMENT="select pipelineQueue.queueID, samples.runID, samples.sampleName, \
samples.coverageID, samples.callerID, assays.assayName, instruments.instrumentName, \
pipelineQueue.status from pipelineQueue \
join samples on samples.sampleID=pipelineQueue.sampleID \
join assays on assays.assayID = samples.assayID \
join instruments on instruments.instrumentID = samples.instrumentID \
where pipelineQueue.status='$DEFAULT_STATUS';"


date_=`date '+%Y-%m-%d %H'`
CURRENTDT=$(echo $date_ | sed -e 's/ /_/g' -e 's/:/_/g' -e 's/-/_/g' )
createSampleFiles "$CURRENTDT" "proton" "gene50"
# createSampleFiles "$CURRENTDT" "proton" "neuro"

INSTRUMENT_ASSAY_PAIR=()

while read -r queueID  runID  sampleName coverageID callerID assay instrument status;
do
	 log "Submitting job to QSUB:
	 QUEUEID - $queueID
	 QUEUE_STATUS - $status
	 RUNID - $runID
	 SAMPLE - $sampleName
	 COVERAGE - $coverageID
	 CALLER - $callerID
	 ASSAY - $assay
	 INSTRUMENT - $instrument
	 ENVIRONMENT - $ENVIRONMENT"

   ####update sampleAnalysisQueue table and set status of this queue to 1 i.e started processing
   updatestatement="update pipelineQueue SET status=1 where queueID = $queueID;"
	 mysql  --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$updatestatement"


	 if [ "$instrument" == "proton" ] && [ "$assay" == "neuro" ] ; then
		 INSTRUMENT_ASSAY_PAIR+=("$instrument/$assay")
		 echo "$queueID;$runID;$sampleName;$coverageID;$callerID;$assay;$instrument;$ENVIRONMENT" >> ${HOME_RUNDIR}proton_neuro_${CURRENTDT}.samples
	 elif [ "$instrument" == "proton" ] && [ "$assay" == "gene50" ] ; then
		 INSTRUMENT_ASSAY_PAIR+=("$instrument/$assay")
		 echo "$queueID;$runID;$sampleName;$coverageID;$callerID;$assay;$instrument;$ENVIRONMENT" >> ${HOME_RUNDIR}proton_gene50_${CURRENTDT}.samples
   fi

done < <(mysql --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$STATUS_QUERY_STATEMENT" -N)

echo "$INSTRUMENT_ASSAY_PAIR"

unique_INSTRUMENT_ASSAY_PAIR=($(printf "%s\n" "${INSTRUMENT_ASSAY_PAIR[@]}" | sort -u))

echo "$unique_INSTRUMENT_ASSAY_PAIR"

for ins_assay_pair in "${unique_INSTRUMENT_ASSAY_PAIR[@]}" ; do

	current_instrument=$(echo $ins_assay_pair | cut -d "/" -f 1)

	current_assay=$(echo $ins_assay_pair | cut -d "/" -f 2)

	qsub -F "-c$CURRENTDT -e$ENVIRONMENT -u$USER -p$PASSWORD" ${HOME_SHELLDIR}submitJob_${current_instrument}_${current_assay}.sh

done


exit

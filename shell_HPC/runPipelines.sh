#!/bin/bash
#===================================================================================
#
# FILE: runPipelines.sh
#
#DESCRIPTION: This script is run by cronjob.It checks the queue table in NGS
#             database and submits a sample qsub based on instrument and assay.
# OPTIONS: see function ’usage’ below
# REQUIREMENTS:
# COMPANY:Houston Methodist Hospital, Molecular Diagnostic Laboratory
# REVISION:
#===================================================================================

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
# functions
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

################################################################################
# initialize Variables
################################################################################

HOME="/home/pipelines/ngs_${ENVIRONMENT}/"
HOME_RUNDIR="${HOME}run_files/"
HOME_SHELLDIR="${HOME}shell/"
DB="ngs_${ENVIRONMENT}"
DATE_=`date '+%Y-%m-%d %H'`
CURRENTDT=$(echo $DATE_ | sed -e 's/ /_/g' -e 's/:/_/g' -e 's/-/_/g' )
DEFAULT_STATUS=0
STATUS_QUERY_STATEMENT="select pipelineQueue.queueID, samples.runID, samples.sampleName, \
samples.coverageID, samples.callerID, assays.assayName, instruments.instrumentName, \
pipelineQueue.status from pipelineQueue \
join samples on samples.sampleID=pipelineQueue.sampleID \
join assays on assays.assayID = samples.assayID \
join instruments on instruments.instrumentID = samples.instrumentID \
where pipelineQueue.status='$DEFAULT_STATUS';"
INSTRUMENT_ASSAY_PAIR=()

################################################################################
#
################################################################################

while read -r queueID runID  sampleName coverageID callerID assay instrument status; do

	 log "Submitting job to QSUB:
	 ENVIRONMENT - $ENVIRONMENT
	 INSTRUMENT - $instrument
	 ASSAY - $assay
	 QUEUEID - $queueID
	 QUEUE_STATUS - $status
	 RUNID - $runID
	 SAMPLE - $sampleName
	 COVERAGE - $coverageID
	 CALLER - $callerID"

   ####update sampleAnalysisQueue table and set status of this queue to 1 i.e started processing
   updatestatement="update pipelineQueue SET status=1 where queueID = $queueID;"
	 mysql  --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$updatestatement"


	 if [ -f ${HOME_RUNDIR}${instrument}_${assay}_${CURRENTDT}.samples ]; then
	 		createSampleFiles "$CURRENTDT" "$instrument" "$assay"
   fi

	 if [ "$instrument" == "proton" ] && [ "$assay" == "neuro" ] ; then
		 INSTRUMENT_ASSAY_PAIR+=("$instrument/$assay")
		 echo "$queueID;$runID;$sampleName;$coverageID;$callerID;$assay;$instrument;$ENVIRONMENT" >> ${HOME_RUNDIR}${instrument}_${assay}_${CURRENTDT}.samples
	 elif [ "$instrument" == "proton" ] && [ "$assay" == "gene50" ] ; then
		 INSTRUMENT_ASSAY_PAIR+=("$instrument/$assay")
		 echo "$queueID;$runID;$sampleName;$coverageID;$callerID;$assay;$instrument;$ENVIRONMENT" >> ${HOME_RUNDIR}${instrument}_${assay}_${CURRENTDT}.samples
   fi

done < <(mysql --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$STATUS_QUERY_STATEMENT" -N)


################################################################################
# submit QSUBS
################################################################################

### get unique instrument - assay pair
unique_INSTRUMENT_ASSAY_PAIR=($(printf "%s\n" "${INSTRUMENT_ASSAY_PAIR[@]}" | sort -u))

for ins_assay_pair in "${unique_INSTRUMENT_ASSAY_PAIR[@]}" ; do

	current_instrument=$(echo $ins_assay_pair | cut -d "/" -f 1)

	current_assay=$(echo $ins_assay_pair | cut -d "/" -f 2)

	qsub -F "-c$CURRENTDT -e$ENVIRONMENT -u$USER -p$PASSWORD" ${HOME_SHELLDIR}submitJob_${current_instrument}_${current_assay}.sh

done

log "Completed cronjob for runPipelines."

exit

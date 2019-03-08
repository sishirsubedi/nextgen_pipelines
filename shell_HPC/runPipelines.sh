#!/bin/bash
#===============================================================================
#
# FILE: runPipelines.sh
#
#DESCRIPTION: This script is run by cronjob.It checks the queue table in NGS
#             database and submits a sample qsub based on instrument and assay.
# OPTIONS: see function display_usuage below
# REQUIREMENTS:
# COMPANY:Houston Methodist Hospital, Molecular Diagnostic Laboratory
# REVISION:
#===============================================================================


# ##############################################################################
# # functions
# ##############################################################################


display_usuage()
{
cat <<EOF >> /dev/stderr

 USAGE: $0

 OPTIONS:
 	-e ENVIRONMENT
 	-u USER
	-p PASSWORD

EOF
}

parse_options()
{
    IMPORT=0

		while getopts "hz:e:u:p:" opt; do

				case $opt in
					h)
					display_usuage
					exit 1
					;;
	        z)
					IMPORT=$OPTARG
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

    if [ $IMPORT -gt 0 ] ; then
        return 0
    fi

    return 1
}

load_modules()
{
      source /home/pipelines/ngs_${ENVIRONMENT}/shell/modules/ngs_utils.sh
}

# ##############################################################################
# workflow:check_db_queue
# ##############################################################################
check_db_queue()
{
	default_Status=0
	status_query_statement="select pipelineQueue.queueID, samples.runID, samples.sampleName, \
	samples.coverageID, samples.callerID, assays.assayName, instruments.instrumentName, \
	pipelineQueue.status from pipelineQueue \
	join samples on samples.sampleID=pipelineQueue.sampleID \
	join assays on assays.assayID = samples.assayID \
	join instruments on instruments.instrumentID = samples.instrumentID \
	where pipelineQueue.status='$default_Status';"

  while read -r queueID runID  sampleName coverageID callerID assay instrument status; do

			log_info "Submitting job to QSUB:
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


		 if [ -f ${HOME_RUN}${instrument}_${assay}_${CURRENTDT}.samples ]; then
		 		create_file "$HOME_RUN"  "${instrument}_${assay}_${CURRENTDT}.samples"
	   fi

		 INSTRUMENT_ASSAY_PAIR+=("$instrument/$assay")
		 echo "$queueID;$runID;$sampleName;$coverageID;$callerID;$assay;$instrument;$ENVIRONMENT" >> ${HOME_RUN}${instrument}_${assay}_${CURRENTDT}.samples


  done < <(mysql --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$status_query_statement" -N)

}

# ##############################################################################
# workflow:submit job
# ##############################################################################
submit_job()
{
### get unique instrument - assay pair
unique_INSTRUMENT_ASSAY_PAIR=($(printf "%s\n" "${INSTRUMENT_ASSAY_PAIR[@]}" | sort -u))

for inst_assay_pair in "${unique_INSTRUMENT_ASSAY_PAIR[@]}" ; do

	current_instrument=$(echo $inst_assay_pair | cut -d "/" -f 1)

	current_assay=$(echo $inst_assay_pair | cut -d "/" -f 2)

  #qsub -F "-c$CURRENTDT -e$ENVIRONMENT -u$USER -p$PASSWORD" ${HOME_SHELL}submitJob_${current_instrument}_${current_assay}.sh

done

log_info "Completed cronjob for $0."

}

# ##############################################################################
# main
# ##############################################################################
main()
{
    parse_options $*

    if [ $? -eq 0 ]
    then
			  echo "Error in previous step. Aborting $0"
        exit 0
    fi

		# ##########################################################################
		# initialize Variables
		# ##########################################################################

    load_modules

		HOME="/home/pipelines/ngs_${ENVIRONMENT}/"
		HOME_RUN="${HOME}run_files/"
		HOME_SHELL="${HOME}shell/"
		DB="ngs_${ENVIRONMENT}"
		DATE_=`date '+%Y-%m-%d %H'`
		CURRENTDT=$(echo $DATE_ | sed -e 's/ /_/g' -e 's/:/_/g' -e 's/-/_/g' )
		INSTRUMENT_ASSAY_PAIR=()

    # ##########################################################################
    # workflows
    # ##########################################################################

    check_db_queue

		submit_job
}


# ##############################################################################
# run main
# ##############################################################################

main $*

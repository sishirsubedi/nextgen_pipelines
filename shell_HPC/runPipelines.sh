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


		 if [ ! -f ${HOME_RUN}${CURRENTDT}_Queued.samples ]; then
		 		create_file "$HOME_RUN"  "${CURRENTDT}_Queued.samples"
        echo "##queueID;runID;sampleName;coverageID;callerID;assay;instrument;ENVIRONMENT##RUNDATE-$CURRENTDT"  >> ${HOME_RUN}${CURRENTDT}_Queued.samples
	   fi
		 echo "$queueID;$runID;$sampleName;$coverageID;$callerID;$assay;$instrument;$ENVIRONMENT" >> ${HOME_RUN}${CURRENTDT}_Queued.samples


  done < <(mysql --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$status_query_statement" -N)

}


# ##############################################################################
# workflow:submit job
# ##############################################################################
submit_job()
{

  tail -n +2 ${HOME_RUN}${CURRENTDT}_Queued.samples | while IFS=';' read -ra line; do

    queueID="${line[0]}"
  	runID="${line[1]}"
  	sampleName="${line[2]}"
  	coverageID="${line[3]}"
  	callerID="${line[4]}"
    assay="${line[5]}"
    instrument="${line[6]}"


    home_analysis_instrument="${HOME_ANALYSIS}${instrument}Analysis/"


    if [ "$instrument" == "proton" ] ; then

      runfolder=$(ls -d /home/${instrument}/*$runID)
      runname=${runfolder##*/}

      create_dir ${home_analysis_instrument}$runname
      create_dir ${home_analysis_instrument}${runname}/$sampleName

      working_dir="${home_analysis_instrument}${runname}/${sampleName}/"

      /opt/torque/bin/qsub -d ${working_dir}  \
           -F "-r$runID -s$sampleName -c$coverageID -v$callerID -a$assay -i$instrument -e$ENVIRONMENT -q$queueID -u$USER -p$PASSWORD" \
           ${HOME_SHELL}ionPipelineInterface.sh

    fi

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
			  log_error "Import flag non-zero. Aborting $0"
        exit 1
    fi

		# ##########################################################################
		# initialize Variables
		# ##########################################################################

    load_modules

		HOME="/home/pipelines/ngs_${ENVIRONMENT}/"
		HOME_RUN="${HOME}run_files/"
		HOME_SHELL="${HOME}shell/"
    HOME_ANALYSIS="/home/environments/ngs_${ENVIRONMENT}/"
		DB="ngs_${ENVIRONMENT}"
		CURRENTDT=`date '+%Y_%m_%d_%H_%M_%S'`


    # ##########################################################################
    # workflows
    # ##########################################################################

    check_db_queue

    if [ ! -f  ${HOME_RUN}${CURRENTDT}_Queued.samples  ] ; then
      exit 0
    fi

		submit_job
}


# ##############################################################################
# run main
# ##############################################################################

main $*

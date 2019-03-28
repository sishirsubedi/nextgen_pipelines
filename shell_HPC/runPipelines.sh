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
      source /home/pipelines/ngs_${ENVIRONMENT}/shell/modules/runBcl2fastqPipeline.sh
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

	   ###update sampleAnalysisQueue table and set status of this queue to 1 i.e started processing
	   updatestatement="update pipelineQueue SET status=1 where queueID = $queueID;"
		 mysql  --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$updatestatement"

     if [ "$instrument" == "proton" ] ; then

       if [ ! -f ${HOME_RUN}${CURRENTDT}_ProtonQueued.samples ]; then
         create_file "$HOME_RUN"  "${CURRENTDT}_ProtonQueued.samples"
         echo "##queueID;runID;sampleName;coverageID;callerID;assay;instrument;ENVIRONMENT##RUNDATE-$CURRENTDT"  >> ${HOME_RUN}${CURRENTDT}_ProtonQueued.samples
       fi
       echo "$queueID;$runID;$sampleName;$coverageID;$callerID;$assay;$instrument;$ENVIRONMENT" >> ${HOME_RUN}${CURRENTDT}_ProtonQueued.samples

     elif [ "$instrument" == "nextseq" ] ;then

       if [ ! -f ${HOME_RUN}${CURRENTDT}_IlluminaQueued.samples ]; then
         create_file "$HOME_RUN"  "${CURRENTDT}_IlluminaQueued.samples"
         echo "##queueID;runID;sampleName;coverageID;callerID;assay;instrument;ENVIRONMENT##RUNDATE-$CURRENTDT"  >> ${HOME_RUN}${CURRENTDT}_IlluminaQueued.samples
       fi

       echo "$queueID;$runID;$sampleName;$coverageID;$callerID;$assay;$instrument;$ENVIRONMENT" >> ${HOME_RUN}${CURRENTDT}_IlluminaQueued.samples

     fi


  done < <(mysql --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$status_query_statement" -N)

}


# ##############################################################################
# workflow:submit job
# ##############################################################################
submit_jobs_proton()
{

  if [ ! -f  ${HOME_RUN}${CURRENTDT}_ProtonQueued.samples  ] ; then
    log_info "No samples queued on proton."
    return 1
  fi

  log_info "Samples queued on proton."

  tail -n +2 ${HOME_RUN}${CURRENTDT}_ProtonQueued.samples | while IFS=';' read -ra line; do

    queueID="${line[0]}"
  	runID="${line[1]}"
  	sampleName="${line[2]}"
  	coverageID="${line[3]}"
  	callerID="${line[4]}"
    assay="${line[5]}"
    instrument="${line[6]}"


    echo $queueID $runID $instrument $assay $sampleName

    home_analysis_instrument="${HOME_ANALYSIS}${instrument}Analysis/"


    runFolder=$(ls -d /home/${instrument}/*$runID)
    runName=${runFolder##*/}

    create_dir ${home_analysis_instrument}$runName
    create_dir ${home_analysis_instrument}${runName}/$sampleName

    working_dir="${home_analysis_instrument}${runName}/${sampleName}/"

    /opt/torque/bin/qsub -d ${working_dir}  \
           -F "-r$runID -s$sampleName -c$coverageID -v$callerID -a$assay -i$instrument -e$ENVIRONMENT -q$queueID -u$USER -p$PASSWORD" \
           ${HOME_SHELL}protonPipelineInterface.sh

  done
}

submit_sample_illumina()
{

  queueID=$1
  runID=$2
  instrument=$3
  assay=$4
  sampleName=$5

  runFolder=$(ls -d /home/$instrument/*_"$runID"_*)
  runName=${runFolder##/home/$instrument/}

  home_analysis_instrument="${HOME_ANALYSIS}${instrument}Analysis/"

  create_dir ${home_analysis_instrument}$runName
  create_dir ${home_analysis_instrument}${runName}/$sampleName

  working_dir="${home_analysis_instrument}${runName}/${sampleName}/"

  log_info "Submitting illuminaPipelineInterface.sh"
  /opt/torque/bin/qsub -d ${working_dir}  \
       -F "-r$runID -s$sampleName -a$assay -i$instrument -e$ENVIRONMENT -q$queueID -u$USER -p$PASSWORD" \
       ${HOME_SHELL}illuminaPipelineInterface.sh
}

submit_jobs_illumina()
{

  if [ ! -f  ${HOME_RUN}${CURRENTDT}_IlluminaQueued.samples  ] ; then
    log_info "No samples queued on illumina."
    return 1
  fi

  log_info "Samples queued on illumina."

  tail -n +2 ${HOME_RUN}${CURRENTDT}_IlluminaQueued.samples | while IFS=';' read -ra line; do

    queueID="${line[0]}"
  	runID="${line[1]}"
  	sampleName="${line[2]}"
  	coverageID="${line[3]}"
  	callerID="${line[4]}"
    assay="${line[5]}"
    instrument="${line[6]}"

    echo $queueID   $runID $instrument $assay $sampleName
    fastqStatus=$(mysql --user="$USER" --password="$PASSWORD" --database="$DB" -se "select status from pipelineStatusBcl2Fastq where runID='$runID'")

    if [ $fastqStatus == "0" ] ; then

      lockdir=${HOME_RUN}nextseq_${runID}.lock

      ## mkdir is atomic, file or file variable is not
      if mkdir "$lockdir" ; then

        log_info "Acquired lock for $lockdir"
        log_info "bcl2fastq_running_now"
        #updateStatus "$queueID_re" "bcl2fastq_running_now" "$environment" "$user"  "$password"

        run_bcl2fastq $USER  $PASSWORD  $DB  $DB_HOST  $runID $ENVIRONMENT

        ### sleep for an hour and keep checking database if bcl2fastq is completed
        checkfastqStatus="0"
        while [ $checkfastqStatus == "0" ]; do

            sleep 1h

            checkfastqStatus=$(mysql --user="$USER" --password="$PASSWORD" --database="$DB" -se "select status from pipelineStatusBcl2Fastq where runID='$runID'")

            log_info "checking bcl2fastq status -- $checkfastqStatus"

        done


        #updateStatus "$queueID_re" "bcl2fastq_completed_now" "$environment" "$user"  "$password"
        submit_sample_illumina $queueID  $runID $instrument $assay $sampleName

        rm -rf "$lockdir"

      else

        tail -n +2 ${HOME_RUN}${CURRENTDT}_IlluminaQueued.samples | while IFS=';' read -ra line; do

          queueID_re="${line[0]}"
          runID_re="${line[1]}"

    		      if [ $runID_re == $runID ]; then

    			           updatestatement="UPDATE pipelineQueue SET pipelineQueue.status=0 WHERE queueID = $queueID_re;"
    			           mysql --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$updatestatement"

    			           #updateStatus "$queueID_re" "bcl2fastq_wait" "$environment" "$user"  "$password"
    		      fi

    	  done
      fi
    else

      #updateStatus "$queueID_re" "bcl2fastq_completed_past" "$environment" "$user"  "$password"
      submit_sample_illumina $queueID  $runID $instrument $assay $sampleName

    fi

  done

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
    DB_HOST="hhplabngsp01"
		CURRENTDT=`date '+%Y_%m_%d_%H_%M_%S'`

    # ##########################################################################
    # workflows
    # ##########################################################################

    check_db_queue

		submit_jobs_proton

    submit_jobs_illumina

    log_info "Completed cronjob for $0."
}


# ##############################################################################
# run main
# ##############################################################################

main $*

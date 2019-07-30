#!/bin/bash
#PBS -l nodes=1
#PBS -l walltime=99:00:00
#PBS -q default
#PBS -k eo

# ##############################################################################
# # functions
# ##############################################################################

display_usage()
{
cat <<EOF >> /dev/stderr

 USAGE: $0

 OPTIONS:
 	-d DIR
 	-n NORMAL
	-t TUMOR
  -v VARIANT CALLING
  -e ENVIRONMENT
  -u USER
  -p PASSWORD
EOF
}

parse_options()
{
    IMPORT=0

		while getopts "hz:d:n:t:v:e:u:p:" opt; do

				case $opt in
					h)
					display_usage
					exit 1
					;;
	        z)
					IMPORT=$OPTARG
					;;
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

    if [ $IMPORT -gt 0 ] ; then
        return 0
    fi

    return 1
}

load_modules()
{
	source /home/pipelines/ngs_${ENVIRONMENT}/shell/modules/ngs_utils.sh
}

run_aligner()
{
  log_info "Running Single Sample Aligner only -- TMB Pipeline!"

   bash ${HOME_SHELLDIR}runAlignment_Interface.sh \
	 	 		-s   $NORMAL \
	 	 		-f   $FASTQ_DIR  \
	 	 		-o   $TMB_AL_OUT
}

run_VCsingle()
{
 log_info "Running Single Sample TMB Pipeline!"

 if [ ! -f ${TMB_AL_OUT}${NORMAL}/Alignment/${NORMAL}.sorted.rmdups.bam ] ; then

	 log_info "${NORMAL}.sorted.rmdups.bam file not found"
   log_info "Running Alignment + UNpaired Variant Caller Pipeline on -- ${TMB_AL_OUT}${NORMAL}/Alignment/${NORMAL}.sorted.rmdups.bam"

   bash ${HOME_SHELLDIR}runAlignment_Interface.sh \
	 -s $NORMAL \
	 -f $FASTQ_DIR  \
	 -o $TMB_AL_OUT


	bash /home/pipelines/ngs_test/shell/runVCaller_interface.sh \
  -n ${TMB_AL_OUT}${NORMAL}/Alignment/${NORMAL}.sorted.rmdups.bam \
  -t "NONE" \
  -v varscan-strelka-mutect \
  -o ${TMB_VC_UP_OUT}

else

  log_info "Alignment Completed:-${NORMAL}.sorted.rmdups.bam file found"
  log_info "Running UNpaired Variant Caller Pipeline on -- ${TMB_AL_OUT}${NORMAL}/Alignment/${NORMAL}.sorted.rmdups.bam"

	bash /home/pipelines/ngs_test/shell/runVCaller_interface.sh \
  -n ${TMB_AL_OUT}${NORMAL}/Alignment/${NORMAL}.sorted.rmdups.bam \
  -t "NONE" \
  -v varscan-strelka-mutect \
  -o ${TMB_VC_UP_OUT}

fi

}

run_VCpaired()
{

  log_info "Running NORMAL/TUMOR Pair TMB Pipeline!"

	/opt/parallel/bin/parallel "bash ${HOME_SHELLDIR}runAlignment_Interface.sh \
	    -s   {}  \
	    -f   $FASTQ_DIR  \
	    -o   $TMB_AL_OUT " echo ::: "$NORMAL" "$TUMOR"

	if [ ! -f ${TMB_AL_OUT}${NORMAL}/Alignment/${NORMAL}.sorted.rmdups.bam ] ; then
		log_info "Error: ${NORMAL}.sorted.rmdups.bam file not found"
		exit
	fi

	if [ ! -f ${TMB_AL_OUT}${TUMOR}/Alignment/${TUMOR}.sorted.rmdups.bam ] ; then
		log_info "Error: ${TUMOR}.sorted.rmdups.bam file not found"
		exit
	fi

	bash /home/pipelines/ngs_test/shell/runVCaller_interface.sh \
	    -n ${TMB_AL_OUT}${NORMAL}/Alignment/${NORMAL}.sorted.rmdups.bam \
	    -t ${TMB_AL_OUT}${TUMOR}/Alignment/${TUMOR}.sorted.rmdups.bam \
	    -v varscan-strelka-mutect \
			-o ${TMB_VC_OUT}

}

run_dbUpdate()
{

    tmb_results="$TMB_VC_OUT${TUMOR}_${NORMAL}/${TUMOR}_${NORMAL}.tmb.result"

    tmb_results_statement="load data local infile '$tmb_results' into table sampleTumorMutationBurden FIELDS TERMINATED BY ',' (sampleID,TMBPair,TMBTotalVariants,TMBScore,TMBGroup)"
    mysql --host="$DB_HOST" --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$tmb_results_statement"

}


get_SeqStats()
{

  /opt/python3/bin/python ${HOME}/python/getExomeSeqStat.py  "${TMB_AL_OUT}/${TUMOR}/Alignment"  "${TMB_AL_OUT}/${NORMAL}/Alignment"

}
# ##############################################################################
# main
# ##############################################################################
main()
{
	parse_options $*

	if [ $? -eq 0 ] ; then
			log_error "Import flag non-zero. Aborting $0"
			exit 1
	fi


	# ##########################################################################
	# initialize Variables
	# ##########################################################################

  load_modules


	FASTQ_DIR="/home/nextseq/${DIR}/out1/"
	TMB_OUT="/home/environments/ngs_${ENVIRONMENT}/exomeAnalysis/tmbAnalysis/${DIR}/"

	create_dir ${TMB_OUT}
	TMB_AL_OUT="${TMB_OUT}Single/"
	TMB_VC_OUT="${TMB_OUT}Paired/"
	TMB_VC_UP_OUT="${TMB_OUT}UNpaired/"

	create_dir ${TMB_AL_OUT}
	create_dir ${TMB_VC_OUT}
	create_dir ${TMB_VC_UP_OUT}


	INSTRUMENT="nextseq"
	ASSAY="TMB"
	HOME="/home/pipelines/ngs_${ENVIRONMENT}/"
	HOME_RUNDIR="${HOME}run_files/"
	HOME_SHELLDIR="${HOME}shell/"
	DB="ngs_${ENVIRONMENT}"
  DB_HOST="hhplabngsp01"

  LOG_FILE="${TMB_OUT}${NORMAL}_${TUMOR}.log"

  exec >  >(tee -a ${LOG_FILE})
  exec 2> >(tee -a ${LOG_FILE} >&2)

  log_info "Starting TMB interface"

  show_pbsinfo

	##############################################################################
	log_info "Parameters are-
Directory-$DIR
Normal Sample-$NORMAL
Tumor Sample-$TUMOR
Variant Calling-$VC
Environment- $ENVIRONMENT
Fastq Directory-$FASTQ_DIR
TMB Analysis Directory-$TMB_OUT
	$TMB_AL_OUT
	$TMB_VC_OUT
	$TMB_VC_UP_OUT
Instrument-$INSTRUMENT
Assay-$ASSAY
Home-$HOME
Home Run Files-$HOME_RUNDIR
Home Shell Files-$HOME_SHELLDIR
Database-$DB"
	##############################################################################


	if [[ $VC == "NO"  &&  $TUMOR == "NO" ]]; then

    log_info "Running Single Sample Aligner only -- TMB Pipeline!"
    run_aligner

  elif [[ $VC == "YES"  &&  $TUMOR == "NO" ]]; then

    log_info "Running Single Sample TMB Pipeline!"
		run_VCsingle

  elif [[ $VC == "YES"  &&  $TUMOR != "NO" ]]; then

    log_info "Running NORMAL/TUMOR Pair TMB Pipeline!"
    # run_VCpaired
    # run_dbUpdate

  else

    log_info "Terminated !! Incorrect parameter for TMB Pipeline!"

  fi

}

# ##############################################################################
# run main
# ##############################################################################
main $*
TMB

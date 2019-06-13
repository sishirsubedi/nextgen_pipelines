#!/bin/bash
#PBS -l nodes=1:ppn=8
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

EOF
}

parse_options()
{
    IMPORT=0

		while getopts "hz:d:n:t:e:" opt; do

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
					e)
					ENVIRONMENT=$OPTARG
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

run_paired()
{
	/opt/parallel/bin/parallel "bash ${HOME_SHELLDIR}runAlignment_Interface.sh \
	    -s   {}  \
	    -f   $FASTQ_DIR  \
	    -o   $EXOME_AL_OUT " echo ::: "$NORMAL" "$TUMOR"

	if [ ! -f ${EXOME_AL_OUT}${NORMAL}/Alignment/${NORMAL}.sorted.rmdups.bam ] ; then
		echo "Error: ${NORMAL}.sorted.rmdups.bam file not found"
		exit
	fi

	if [ ! -f ${EXOME_AL_OUT}${TUMOR}/Alignment/${TUMOR}.sorted.rmdups.bam ] ; then
		echo "Error: ${TUMOR}.sorted.rmdups.bam file not found"
		exit
	fi

	bash /home/pipelines/ngs_test/shell/runVCaller_interface.sh \
	    -n ${EXOME_AL_OUT}${NORMAL}/Alignment/${NORMAL}.sorted.rmdups.bam \
	    -t ${EXOME_AL_OUT}${TUMOR}/Alignment/${TUMOR}.sorted.rmdups.bam \
	    -v varscan-strelka-mutect \
			-o ${EXOME_VC_OUT}
}

run_single()
{
 echo "Running Single Sample Exome Pipeline!"

 if [ ! -f ${EXOME_AL_OUT}${NORMAL}/Alignment/${NORMAL}.sorted.rmdups.bam ] ; then
	 echo "Error: ${NORMAL}.sorted.rmdups.bam file not found"
	 bash ${HOME_SHELLDIR}runAlignment_Interface.sh \
	 	 		-s   $NORMAL \
	 	 		-f   $FASTQ_DIR  \
	 	 		-o   $EXOME_AL_OUT
else
	echo "Running UNpaired VC on -- ${EXOME_AL_OUT}${NORMAL}/Alignment/${NORMAL}.sorted.rmdups.bam"

	bash /home/pipelines/ngs_test/shell/runVCallerUNpaired_interface.sh \
	 	 -s ${EXOME_AL_OUT}${NORMAL}/Alignment/${NORMAL}.sorted.rmdups.bam \
	 	 -v varscan-strelka-mutect \
	 	 -o ${EXOME_VC_UP_OUT}
fi

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
	EXOME_OUT="/home/environments/ngs_${ENVIRONMENT}/exomeAnalysis/${DIR}/"

	create_dir ${EXOME_OUT}
	EXOME_AL_OUT="${EXOME_OUT}Single/"
	EXOME_VC_OUT="${EXOME_OUT}Paired/"
	EXOME_VC_UP_OUT="${EXOME_OUT}UNpaired/"

	create_dir ${EXOME_AL_OUT}
	create_dir ${EXOME_VC_OUT}
	create_dir ${EXOME_VC_UP_OUT}


	INSTRUMENT="nextseq"
	ASSAY="exome"
	HOME="/home/pipelines/ngs_${ENVIRONMENT}/"
	HOME_RUNDIR="${HOME}run_files/"
	HOME_SHELLDIR="${HOME}shell/"
	DB="ngs_${ENVIRONMENT}"

	##############################################################################
	echo "Parameters are-
	$DIR
	$NORMAL
	$TUMOR
	$ENVIRONMENT
	$FASTQ_DIR
	$EXOME_OUT
	$EXOME_AL_OUT
	$EXOME_VC_OUT
	$EXOME_VC_UP_OUT
	$INSTRUMENT
	$ASSAY
	$HOME
	$HOME_RUNDIR
	$HOME_SHELLDIR
	$DB"
	##############################################################################


	if [[ $TUMOR == "NONE" ]]; then
  	echo "Running Single Sample Exome Pipeline!"
		run_single
	else
		run_paired
  fi

}

# ##############################################################################
# run main
# ##############################################################################
main $*

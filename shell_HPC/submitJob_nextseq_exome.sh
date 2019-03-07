##############################################################################
#
# Houston Methodist Hospital
# Molecular Diagnostic
#
#Description:
#This script checks samples queued in instrument/assay specific txt file
# and calls appropriate interface script.
#Allocates appropriate PBS parameters.
##############################################################################

#!/bin/bash
#PBS -l nodes=n002.cluster.com
#PBS -l walltime=99:00:00
#PBS -q default
#PBS -k eo

################################################################################
# assign Variables
################################################################################

# while getopts :c:e:u:p: opt; do
# 				case $opt in
#         c)
#   			FILE_ID=$OPTARG
#   			;;
# 			  e)
# 				ENVIRONMENT=$OPTARG
# 				;;
# 			  u)
# 			  USER=$OPTARG
# 			  ;;
# 			  p)
# 			  PASSWORD=$OPTARG
# 				;;
# 			  :)
# 				echo "Option -$OPTARG requires an argument."
# 				;;
# 			  \?)
# 				echo "Invalid option: -$OPTARG"
# 		esac
# done
# 	  shift $((OPTIND -1))


################################################################################
# functions
################################################################################
function log() {
	MESSAGE=$1
	TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
	SCRIPT="submitJob_nextseq_exome"
	echo " [ $TIMESTAMP ] [ $SCRIPT ] : $MESSAGE "
}


################################################################################

FASTQ_DIR="/home/nextseq/190128_NS500761_0287_AHLW2JBGX9/out1/"
EXOME_AL_OUT="/home/environments/ngs_test/exomeAnalysis/Single/"
EXOME_VC_OUT="/home/environments/ngs_test/exomeAnalysis/Paired/"
NORMAL="COLO-829BL_Win10_S4"
TUMOR="COLO-829_Win10_S5"
ENVIRONMENT="test"
################################################################################

################################################################################
# initialize Variables
################################################################################
INSTRUMENT="nextseq"
ASSAY="exome"
HOME="/home/pipelines/ngs_${ENVIRONMENT}/"
HOME_RUNDIR="${HOME}run_files/"
HOME_SHELLDIR="${HOME}shell/"
DB="ngs_${ENVIRONMENT}"

################################################################################
#
################################################################################

/opt/parallel/bin/parallel "bash ${HOME_SHELLDIR}runAlignment_Interface_window.sh \
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
    -v varscan \
    -o ${EXOME_VC_OUT}

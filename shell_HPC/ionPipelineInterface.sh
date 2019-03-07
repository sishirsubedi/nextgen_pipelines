##############################################################################
#
# Houston Methodist Hospital
# Molecular Diagnostic
#
#Description:
#This script checks parameters, creates run and sample directory
# and calls pipeline script.
##############################################################################

#!/bin/bash

# ##############################################################################
# # functions
# ##############################################################################

display_usuage()
{
cat <<EOF >> /dev/stderr

 USAGE: $0

 OPTIONS:
 r - runID
 s - sampleName
 c - coverageID
 v - callerID
 a - assay
 i - instrument
 e - environment
 q - queueID
 u - user
 p - password

EOF
}

parse_options()
{
    IMPORT=0

		while getopts "hz:r:s:c:v:a:i:e:q:u:p:" opt ; do
				case $opt in
					h)
						 display_usuage
						 exit 1
						 ;;
					z)
						IMPORT=$OPTARG
						;;
					r)
						RUNID=$OPTARG
						;;
					s)
					  SAMPLENAME=$OPTARG
						;;
					c)
					  COVERAGEID=$OPTARG
					  ;;
					v)
						CALLERID=$OPTARG
						;;
					a)
						ASSAY=$OPTARG
						;;
					i)
				    INSTRUMENT=$OPTARG
					  ;;
				  e)
				    ENVIRONMENT=$OPTARG
					  ;;
				  q)
						QUEUEID=$OPTARG
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

log()
{
 MESSAGE=$1
 TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
 SCRIPT=$( basename $0 )
 echo " [ $TIMESTAMP ] [ $SCRIPT ] : $MESSAGE "
}


create_dir()
{
	n_dir=$1
	if [ ! -d $n_dir ] ; then
		mkdir $n_dir
	fi
	chmod 775 $n_dir
}

run_gene50()
{
	NEURO_EXCLUDED_AMPLICON="/home/doc/ref/neuralRef/excludedAmplicon.txt"
	NEURO_EXCLUDED_DESIGN="/home/doc/ref/neuralRef/IAD87786_179_Designed.excluded.bed"

	bash ${HOME_SHELLDIR}ionPipeline.sh -r $RUNID -s $SAMPLENAME -c $COVERAGEID -v $CALLERID \
																		-i $INSTRUMENT  -e $NEURO_EXCLUDED_AMPLICON -a $NEURO_EXCLUDED_DESIGN \
																		-n $ENVIRONMENT -q $QUEUEID -u $USER -p $PASSWORD
}

run_neuro()
{
	bash ${HOME_SHELLDIR}ionPipeline.sh -r $RUNID -s $SAMPLENAME -c $COVERAGEID -v $CALLERID \
																		-i $INSTRUMENT \
																		-n $ENVIRONMENT -q $QUEUEID -u $USER -p $PASSWORD
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

		################################################################################
		# initialize variables
		################################################################################

		VARIANTFOLDER=$(ls -d /home/$INSTRUMENT/*$runID/plugin_out/"$CALLERID")
		AMPLICONFOLDER=$(ls -d /home/$INSTRUMENT/*"$RUNID"/plugin_out/"$COVERAGEID")
		RUNFOLDER=$(ls -d /home/$INSTRUMENT/*$RUNID)
		RUNNAME=${RUNFOLDER##*/}
		HOME="/home/environments/ngs_${ENVIRONMENT}/${INSTRUMENT}Analysis/"
		HOME_SHELLDIR="/home/pipelines/ngs_${ENVIRONMENT}/shell/"


		log " Running ionPipeline Interface for :
		assay : $ASSAY
		instrument : $INSTRUMENT
		runID : $RUNID
		sampleName : $SAMPLENAME
		coverageID : $COVERAGEID
		callerID : $CALLERID
		environment : $ENVIRONMENT
		queueID : $QUEUEID "

		create_dir ${HOME}$RUNNAME
		create_dir ${HOME}${RUNNAME}/$SAMPLENAME
		create_dir ${HOME}${RUNNAME}/${SAMPLENAME}/$CALLERID
		create_dir ${HOME}${RUNNAME}/${SAMPLENAME}/$COVERAGEID

		exec >  >(tee -a ${HOME}${RUNNAME}/${SAMPLENAME}/process.log)
		exec 2> >(tee -a ${HOME}${RUNNAME}/${SAMPLENAME}/process.log >&2)

		if [ $assay == "neuro" ] ; then
			run_gene50
		elif [ $assay == "gene50" ] ; then
			run_neuro
		fi
}


# ##############################################################################
# run main
# ##############################################################################

main $*

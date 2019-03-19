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
#PBS -l nodes=1
#PBS -l walltime=2:00:00
#PBS -q default
#PBS -o ${PBS_JOBNAME}.out
#PBS -e ${PBS_JOBNAME}.err
#PBS -k eo

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

load_modules()
{
      source /home/pipelines/ngs_${ENVIRONMENT}/shell/modules/ngs_utils.sh
}


create_rundate()
{
  declare -A months
  months=( ["Jan"]="01" ["Feb"]="02" ["Mar"]="03" ["Apr"]="04" ["May"]="05" ["Jun"]="06" ["Jul"]="07" ["Aug"]="08" ["Sep"]="09" ["Oct"]="10" ["Nov"]="11" ["Dec"]="12" )
  if [ -f /home/${INSTRUMENT}/*"$RUNID"/InitLog.txt ]
  then
  	runDate=$(head -n 1 /home/${INSTRUMENT}/*"$RUNID"/InitLog.txt)
  	year1=$(echo $runDate |cut -d ' ' -f 5)
  	year=${year1%:}
  	day=$(echo $runDate |cut -d ' ' -f 3)
  	monthWord=$(echo $runDate |cut -d ' ' -f 2)
  	month=${months["$monthWord"]}
  	date=$year-$month-$day
  	echo $date > ${HOME}${RUNNAME}/${SAMPLENAME}/runDate.txt
  else
  	log_info "Warning: InitLog.txt SAMPLE_FILE not found, run date will not be entered"
  fi
}


prep_neuro()
{
  NEURO_EXCLUDED_DESIGN="/home/doc/ref/neuralRef/IAD87786_179_Designed.excluded.bed"

  #bedtools -u flag to write original A entry once if any overlaps found in B
  /opt/software/bedtools-2.17.0/bin/bedtools intersect -u -a $PROTON_VCF -b $NEURO_EXCLUDED_DESIGN > ${HOME}${RUNNAME}/${SAMPLENAME}/${CALLERID}/TSVC_variants.filter.vcf


  NEURO_EXCLUDED_AMPLICON="/home/doc/ref/neuralRef/excludedAmplicon.txt"
  # -v means "invert the match" in grep, in other words, return all non matching lines.
  grep -v -f $NEURO_EXCLUDED_AMPLICON $PROTON_AMPLICON > ${HOME}${RUNNAME}/${SAMPLENAME}/${COVERAGEID}/amplicon.filter.txt

}

prep_gene50()
{

  ln -s $PROTON_VCF ${HOME}${RUNNAME}/${SAMPLENAME}/${CALLERID}/TSVC_variants.filter.vcf

  ln -s $PROTON_AMPLICON ${HOME}${RUNNAME}/${SAMPLENAME}/${COVERAGEID}/amplicon.filter.txt

}

run_ionPipeline()
{

  update_status "$QUEUEID" "Started" "$ENVIRONMENT" "$USER"  "$PASSWORD"

  bash ${HOME_SHELLDIR}ionPipeline.sh -d "${HOME}${RUNNAME}/${SAMPLENAME}/" \
        -s $SAMPLENAME -c $COVERAGEID -v $CALLERID -e $ENVIRONMENT -q $QUEUEID -u $USER -p $PASSWORD

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

		############################################################################
		# initialize variables
		############################################################################

    load_modules

		VARIANTFOLDER=$(ls -d /home/${INSTRUMENT}/*${RUNID}/plugin_out/"$CALLERID")
    PROTON_VCF="$VARIANTFOLDER/${SAMPLENAME}/TSVC_variants.vcf"

    AMPLICONFOLDER=$(ls -d /home/${INSTRUMENT}/*${RUNID}/plugin_out/"$COVERAGEID")
    PROTON_AMPLICON=$(ls $AMPLICONFOLDER/${SAMPLENAME}/*.amplicon.cov.xls)

    RUNFOLDER=$(ls -d /home/${INSTRUMENT}/*$RUNID)
		RUNNAME=${RUNFOLDER##*/}

    HOME="/home/environments/ngs_${ENVIRONMENT}/${INSTRUMENT}Analysis/"
		HOME_SHELLDIR="/home/pipelines/ngs_${ENVIRONMENT}/shell/"


		create_dir ${HOME}${RUNNAME}/${SAMPLENAME}/$CALLERID
		create_dir ${HOME}${RUNNAME}/${SAMPLENAME}/$COVERAGEID

		exec >  >(tee -a ${HOME}${RUNNAME}/${SAMPLENAME}/process.log)
		exec 2> >(tee -a ${HOME}${RUNNAME}/${SAMPLENAME}/process.log >&2)

    show_pbsinfo

    log_info " Running ionPipeline Interface BY qsub for :
		assay : $ASSAY
		instrument : $INSTRUMENT
		runID : $RUNID
		sampleName : $SAMPLENAME
		coverageID : $COVERAGEID
		callerID : $CALLERID
		environment : $ENVIRONMENT
		queueID : $QUEUEID "




    create_rundate

		if [ $ASSAY == "neuro" ] ; then
			prep_neuro
		elif [ $ASSAY == "gene50" ] ; then
			prep_gene50
		fi

    run_ionPipeline
}


# ##############################################################################
# run main
# ##############################################################################


main $*

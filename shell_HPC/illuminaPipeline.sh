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
 d - sample directory under analysis environment
 s - sampleName
 c - coverageID
 v - callerID
 e - environment
 q - queueID
 u - user
 p - password

EOF
}


parse_options()
{
    IMPORT=0

		while getopts "hz:d:s:e:q:u:p:" opt ; do
				case $opt in
					h)
						 display_usuage
						 exit 1
						 ;;
					z)
						IMPORT=$OPTARG
						;;
					d)
						HOME=$OPTARG
						;;
					s)
						SAMPLENAME=$OPTARG
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
      source /home/pipelines/ngs_${ENVIRONMENT}/shell/modules/ngs_vep.sh
}



run_vep()
{
	log_info "Running VEP"

  start_vep "${HOME}variantAnalysis/${SAMPLENAME}.filter.vcf"  "${HOME}variantAnalysis/${SAMPLENAME}.filter.vep.vcf"

  log_info "Completed VEP"
}

parse_vep()
{
  log_info "Parse VEP"

	python ${HOME_PYTHON}parseVEP.py \
					parseIllumina  \
					-I ${HOME}variantAnalysis/${SAMPLENAME}.filter.vep.vcf\
					-o ${HOME}variantAnalysis/${SAMPLENAME}.filter.vep.parse.vcf
}

filter_vep()
{
	log_info "filter VEP"
	shopt -s nocasematch
	if [[ $SAMPLENAME =~ horizon ]]
	then
		awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 1 && $10 != "null" && $11 >=100 && $11 != "null") print}' ${HOME}variantAnalysis/${SAMPLENAME}.filter.vep.parse.vcf \
		> ${HOME}variantAnalysis/${SAMPLENAME}.filter.vep.parse.filter2.vcf
	else
		awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 10 && $10 != "null" && $11 >=100 && $11 != "null") print}' ${HOME}variantAnalysis/${SAMPLENAME}.filter.vep.parse.vcf \
		> ${HOME}variantAnalysis/${SAMPLENAME}.filter.vep.parse.filter2.vcf
	fi

}

filter_amplicon()
{


  log_info "Filter amplicon"

	sampleCol=$(head -n 1 ${HOME}variantAnalysis/${SAMPLENAME}.amplicon.filter.txt | awk -F"\t" -v sample=$SAMPLENAME '{for (i=1; i<= NF; i++) if($i==sample) print i}')

	if [ $sampleCol = "" ]
		then
			echo "amplicon info not found in amplicon file for sample:$SAMPLENAME"
	fi

	cut -f 1,$sampleCol ${HOME}variantAnalysis/${SAMPLENAME}.amplicon.filter.txt > ${HOME}variantAnalysis/${SAMPLENAME}.amplicon.filter.filter2.txt
	awk '$2 < 100' ${HOME}variantAnalysis/${SAMPLENAME}.amplicon.filter.filter2.txt > ${HOME}variantAnalysis/${SAMPLENAME}.amplicon.filter.filter2.lessThan100.txt

}


update_db()
{

  log_info "Updating Database"

  bash ${HOME_SHELL}runDBUpdate.sh -d $HOME -s $SAMPLENAME -c $COVERAGEID -v $CALLERID -q $QUEUEID -e $ENVIRONMENT -u $USER -p $PASSWORD

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

		HOME_PYTHON="/home/pipelines/ngs_${ENVIRONMENT}/python/"
		HOME_SHELL="/home/pipelines/ngs_${ENVIRONMENT}/shell/"
		DB="ngs_${ENVIRONMENT}"

		log_info " Running illumina pipeline for :
		dir : $HOME
		sampleName : $SAMPLENAME
		environment : $ENVIRONMENT
		queueID : $QUEUEID "


    # update_status "$QUEUEID" "RunningVEP" "$ENVIRONMENT" "$USER"  "$PASSWORD"
		#run_vep

    # update_status "$QUEUEID" "CompletedVEP" "$ENVIRONMENT" "$USER"  "$PASSWORD"
		#
		#parse_vep
		#
		#filter_vep
		#
		#filter_amplicon
		#
    # update_status "$QUEUEID" "UpdatingDatabase" "$ENVIRONMENT" "$USER"  "$PASSWORD"
		#
		update_db

}


# ##############################################################################
# run main
# ##############################################################################

main $*

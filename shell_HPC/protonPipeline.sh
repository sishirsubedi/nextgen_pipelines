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

		while getopts "hz:d:s:c:v:e:q:u:p:" opt ; do
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
					c)
						COVERAGEID=$OPTARG
						;;
					v)
						CALLERID=$OPTARG
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


split_protonVCF()
{

	python ${HOME_PYTHON}splitVcf.py \
				-I ${HOME}${CALLERID}/TSVC_variants.filter.vcf \
				-f FAO,FDP,AF \
				-o ${HOME}${CALLERID}/TSVC_variants.filter.split.vcf
}

run_vep()
{
	log_info "Running VEP"

  start_vep "${HOME}${CALLERID}/TSVC_variants.filter.split.vcf"  "${HOME}${CALLERID}/TSVC_variants.filter.split.vep.vcf"

  log_info "Completed VEP"
}

parse_vep()
{
  log_info "Parse VEP"

	python ${HOME_PYTHON}parseVEP.py \
					parseIonNewVarView \
					-I ${HOME}${CALLERID}/TSVC_variants.filter.split.vep.vcf \
					-o ${HOME}${CALLERID}/TSVC_variants.filter.split.vep.parse.txt
}

filter_vep()
{
	log_info "filter VEP"
	shopt -s nocasematch
	if [[ $SAMPLENAME =~ horizon ]]
	then
		awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 1 && $10 != "null" && $11 >=100 && $11 != "null") print}' ${HOME}${CALLERID}/TSVC_variants.filter.split.vep.parse.txt \
		> ${HOME}${CALLERID}/TSVC_variants.filter.split.vep.parse.filter2.txt
	else
		awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 10 && $10 != "null" && $11 >=100 && $11 != "null") print}' ${HOME}${CALLERID}/TSVC_variants.filter.split.vep.parse.txt \
		> ${HOME}${CALLERID}/TSVC_variants.filter.split.vep.parse.filter2.txt
	fi

}

filter_amplicon()
{
  log_info "Filter amplicon"

	awk -F'\t' -v OFS='\t' '{if($10 < 100) print $4,$10}' ${HOME}${COVERAGEID}/amplicon.filter.txt \
	> ${HOME}${COVERAGEID}/amplicon.lessThan100.txt
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

		log_info " Running proton pipeline for :
		director : $HOME
		sampleName : $SAMPLENAME
		coverageID : $COVERAGEID
		callerID : $CALLERID
		environment : $ENVIRONMENT
		queueID : $QUEUEID "

		split_protonVCF

    update_status "$QUEUEID" "RunningVEP" "$ENVIRONMENT" "$USER"  "$PASSWORD"
		run_vep

    update_status "$QUEUEID" "CompletedVEP" "$ENVIRONMENT" "$USER"  "$PASSWORD"

		parse_vep

		filter_vep

		filter_amplicon

    update_status "$QUEUEID" "UpdatingDatabase" "$ENVIRONMENT" "$USER"  "$PASSWORD"

		update_db

}


# ##############################################################################
# run main
# ##############################################################################

main $*
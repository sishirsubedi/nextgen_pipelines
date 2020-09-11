#!/bin/bash
#===============================================================================
#
# FILE: heme_annotationPipeline.sh
#
#DESCRIPTION: This script is run to annonate and filter amplicon and variant files.
# REQUIREMENTS: illuminaPipelineInterface.sh
# COMPANY:Houston Methodist Hospital, Molecular Diagnostic Laboratory
#===============================================================================

#!/bin/bash

# ##############################################################################
# # functions
# ##############################################################################
display_usage()
{
cat <<EOF >> /dev/stderr

 USAGE: $0

 OPTIONS:
 d - sample directory under analysis environment
 s - sampleName
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
						 display_usage
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

run_snpeff_vep()
{
	log_info "Running SNPEFF"

  snpeff_germline  "${HOME_ANALYSIS}variantAnalysis/${SAMPLENAME}.mutect.combine.v4.vcf"  "${HOME_ANALYSIS}variantAnalysis/${SAMPLENAME}.mutect.snpeff.vcf"

  vep_germline  "${HOME_ANALYSIS}variantAnalysis/${SAMPLENAME}.mutect.combine.v4.vcf"  "${HOME_ANALYSIS}variantAnalysis/${SAMPLENAME}.mutect.vep.vcf"

  log_info "Completed SNPEFF"
}

parse_snpeff_vep_add_annotation()
{
  log_info "Parse SNPEFF"

  /opt/python3/bin/python3 ${HOME_PYTHON}germline_parse_snpeff.py \
  ${HOME}variantAnalysis/${SAMPLENAME}.mutect.snpeff.vcf \
  ${HOME}variantAnalysis/${SAMPLENAME}.mutect.snpeff.parse.vcf \
  ${HOME}variantAnalysis/${SAMPLENAME}.mutect.snpeff.parse_all.vcf

  log_info "Parse VEP"

  /opt/python3/bin/python3 ${HOME_PYTHON}germline_parse_vep.py \
  ${HOME}variantAnalysis/${SAMPLENAME}.mutect.vep.vcf \
  ${HOME}variantAnalysis/${SAMPLENAME}.mutect.vep.parse.vcf


  ### add_annotation
  /opt/python3/bin/python3 ${HOME_PYTHON}germline_add_annotation.py \
  ${SAMPLENAME} \
  ${HOME}variantAnalysis/ \
  ${ENVIRONMENT} ${DB_HOST} ${DB} ${USER}  ${PASSWORD}

}

filter_amplicon()
{

  log_info "Filter amplicon"

 sampleCol=$(head -n 1  ${HOME}variantAnalysis/${SAMPLENAME}.amplicon.filter.txt | awk -F"\t" '{print NF; exit}')

	if [ $sampleCol != "2" ]
		then
			echo "amplicon info not found in amplicon file for sample:$SAMPLENAME"
	fi

	cut -f 1,$sampleCol ${HOME}variantAnalysis/${SAMPLENAME}.amplicon.filter.txt > ${HOME}variantAnalysis/${SAMPLENAME}.amplicon.filter.filter2.txt
	awk '$2 < 100' ${HOME}variantAnalysis/${SAMPLENAME}.amplicon.filter.filter2.txt > ${HOME}variantAnalysis/${SAMPLENAME}.amplicon.filter.filter2.lessThan100.txt

}

update_db()
{

  log_info "Updating Database"

  bash ${HOME_SHELL}heme_proton_common_DBUpdate.sh -d $HOME -s $SAMPLENAME -c "-" -v "-" -q $QUEUEID -e $ENVIRONMENT -u $USER -p $PASSWORD

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
    DB_HOST="hhplabngsp01"

		log_info " Running illumina pipeline for :
		dir : $HOME
		sampleName : $SAMPLENAME
		environment : $ENVIRONMENT
		queueID : $QUEUEID "


    update_status "$QUEUEID" "RunningVEP" "$DB" "$USER"  "$PASSWORD"

		# run_snpeff_vep

    update_status "$QUEUEID" "CompletedVEP" "$DB" "$USER"  "$PASSWORD"

		parse_snpeff_vep_add_annotation

    update_status "$QUEUEID" "UpdatingDatabase" "$DB" "$USER"  "$PASSWORD"

		update_db

}

# ##############################################################################
# run main
# ##############################################################################
main $*

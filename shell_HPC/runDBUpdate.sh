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
 q - queueID
 e - environment
 u - user
 p - password

EOF
}

parse_options()
{
    IMPORT=0

		while getopts "hz:d:s:c:v:q:e:u:p:" opt ; do
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
          q)
  					QUEUEID=$OPTARG
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

get_IDs()
{
  getids_statement="select samples.sampleID, samples.runID, instruments.instrumentName from pipelineQueue \
  join samples on samples.sampleID=pipelineQueue.sampleID \
  join instruments on instruments.instrumentID = samples.instrumentID \
  where pipelineQueue.queueID='$QUEUEID';"

  while  read -r sampleID runID instrument ; do

     SAMPLEID="$sampleID"
     RUNID="$runID"
     INSTRUMENT="$instrument"

  done < <(mysql --host="$DB_HOST" --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$getids_statement" -N)

}

load_variants()
{
    if [ "$INSTRUMENT" == "proton" ] ; then

       variantfile=$(ls ${HOME}${CALLERID}/TSVC_variants.filter.split.vep.parse.filter2.txt)

  		 if [ ! -f $variantfile ]; then
  		 		  update_status "$QUEUEID" "ERROR:VariantFile" "$ENVIRONMENT" "$USER"  "$PASSWORD"
            log_error "Analsysis variant file not found !"
  					exit 1
  		 fi

  		 variantstatement="load data local infile '$variantfile' into table sampleVariants (gene, exon, chr, pos, ref, alt, impact, type, quality, altFreq, readDepth, altReadDepth, consequence, Sift, PolyPhen, HGVSc, HGVSp, dbSNPID, pubmed) set sampleID = '$SAMPLEID'"

       mysql --host="$DB_HOST" --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$variantstatement"
    fi
}

load_amplicons()
{
    /opt/python3/bin/python3 ${HOME_PYTHON}loadSampleAmplicons.py \
    -i ${HOME}${COVERAGEID}/amplicon.filter.txt \
    -o ${HOME}${COVERAGEID}/amplicon.filter.v2.txt \
    -s $SAMPLEID \
    -n $INSTRUMENT

    ampliconfile=$(ls ${HOME}${COVERAGEID}/amplicon.filter.v2.txt)

  	if [ ! -f $ampliconfile ]; then
           update_status "$QUEUEID" "ERROR:AmpliconFile" "$ENVIRONMENT" "$USER"  "$PASSWORD"
  		 		 log_error "Analsysis amplicon file not found !"
  		 		exit
  	fi

   ampliconstatement="load data local infile '$ampliconfile' into table sampleAmplicons FIELDS TERMINATED BY '\t' (sampleID, gene, ampliconName,readDepth )"
   mysql --host="$DB_HOST" --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$ampliconstatement"

}

log_run()
{
/opt/python3/bin/python3 ${HOME_PYTHON}logPipelineRun.py -q "$QUEUEID" -u "$USER" -p "$PASSWORD" -d "$ENVIRONMENT" -b "$DB_HOST"

}


email_user()
{
  #### get pipeline entering user
  usergetstatement="select email from users where userID in (select enteredBy from samples as t1  join pipelineQueue as t2 on  t1.sampleID=t2.sampleID where t2.queueID ='$QUEUEID');"
  useremail=$(mysql --host="$DB_HOST" --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$usergetstatement" -N)

  log_info "Emailing $useremail with completed message"

  message=$'Pipeline Completed for-
  				instrument : '$INSTRUMENT' ,
  				runID : '$RUNID' ,
  				sampleName : '$SAMPLENAME' ,
  				coverageID : '$COVERAGEID' ,
  				callerID : '$CALLERID' ,
  				environment : '$ENVIRONMENT' ,
  				queueID : '$QUEUEID' '
  echo $message

  echo $message |  mailx  -s "Pipeline completed for sample:$SAMPLENAME , run:$RUNID" "$useremail"

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

		log_info " Running db update for :
		environment : $ENVIRONMENT
		queueID : $QUEUEID
    directory : $HOME "

    SAMPLEID=""
    RUNID=""
    INSTRUMENT=""

    get_IDs

    load_variants

    load_amplicons

    log_run

    email_user

    update_status "$QUEUEID" "PipelineCompleted" "$ENVIRONMENT" "$USER"  "$PASSWORD"

}


# ##############################################################################
# run main
# ##############################################################################

main $*

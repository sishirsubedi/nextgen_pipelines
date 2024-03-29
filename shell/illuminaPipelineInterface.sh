#!/bin/bash
#===============================================================================
#
# FILE: illuminaPipelineInterface.sh
#
#DESCRIPTION: This script is run to generate required files
#             for heme and tumor mutation burden assays.
# REQUIREMENTS: runPipelines.sh
# COMPANY:Houston Methodist Hospital, Molecular Diagnostic Laboratory
#===============================================================================

#PBS -q default
#PBS Error_Path=${PBS_JOBNAME}.err
#PBS Output_Path=${PBS_JOBNAME}.out
#PBS -k eo

# ##############################################################################
# # functions
# ##############################################################################
display_usage()
{
cat <<EOF >> /dev/stderr

 USAGE: $0

 OPTIONS:
 r - runID
 s - sampleName
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

		while getopts "hz:r:s:a:i:e:q:u:p:" opt ; do
				case $opt in
					h)
						 display_usage
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
      source /storage/apps/pipelines/ngs_${ENVIRONMENT}/shell/modules/ngs_utils.sh
      source /storage/apps/pipelines/ngs_${ENVIRONMENT}/shell/modules/ngs_align.sh
      source /storage/apps/pipelines/ngs_${ENVIRONMENT}/shell/modules/ngs_varscan.sh
}

create_rundate()
{
	##get runDate information##
	dateString=${RUNNAME%%_*}
  echo $dateString
	year=20${dateString:0:2}
	month=${dateString:2:2}
	day=${dateString:4:2}
	dateString=$year-$month-$day

	echo $dateString > ${HOME_ANALYSIS}runDate.txt

}

heme_generate_variantFile()
{

  bash /storage/apps/pipelines/ngs_${ENVIRONMENT}/shell/heme_interface.sh \
      -r $RUNID \
      -s $SAMPLENAME \
      -a $ASSAY \
      -i $INSTRUMENT \
      -e $ENVIRONMENT \
      -q $QUEUEID  \
      -u $USER \
      -p $PASSWORD

  VARIANT_VCF=$(ls ${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.combine.vcf )
}

heme_generate_ampliconFile()
{
      HEME_EXCLUDED_DESIGN="/storage/database/ngs_doc/heme/trusight-myeloid-amplicon-track.excluded.bed"

      /storage/apps/opt/samtools/bin/samtools bedcov  $HEME_EXCLUDED_DESIGN  ${HOME_ANALYSIS}variantCaller/${SAMPLENAME}.sort.bam | awk ' {print $4,"\t",int($13/($8-$7))} ' > ${HOME_ANALYSIS}variantAnalysis/${SAMPLENAME}.samtools.coverageDepth

      AMPLICON_FILE=$(ls ${HOME_ANALYSIS}variantAnalysis/${SAMPLENAME}.samtools.coverageDepth)

      if [ ! -f "${HOME_ANALYSIS}variantAnalysis/${SAMPLENAME}.samtools.coverageDepth" ] ; then


    log_error "ERROR:AmpliconFile output file not found"
    update_status "$QUEUEID" "ERROR:AmpliconFile" "$DB" "$USER"  "$PASSWORD" "$DB_HOST"
    exit

  else
    log_info "Completed AmpliconFile"

  fi

}

heme_pre_annotation()
{

  HEME_EXCLUDED_DESIGN="/storage/database/ngs_doc/heme/trusight-myeloid-amplicon-track.excluded.bed"
  #bedtools -u flag to write original A entry once if any overlaps found in B

  /storage/apps/opt/bedtools/bedtools2_17/bin/bedtools intersect -u -a $VARIANT_VCF -b $HEME_EXCLUDED_DESIGN > ${HOME_ANALYSIS}variantAnalysis/${SAMPLENAME}.filter.vcf


  HEME_EXCLUDED_AMPLICON="/storage/database/ngs_doc/heme/excludedAmplicons.txt"
  # -v means "invert the match" in grep, in other words, return all non matching lines.
  grep -v -f $HEME_EXCLUDED_AMPLICON $AMPLICON_FILE > ${HOME_ANALYSIS}variantAnalysis/${SAMPLENAME}.amplicon.filter.txt

}

heme_run_annotation()
{
  bash ${HOME_SHELLDIR}heme_annotationPipeline.sh -d $HOME_ANALYSIS \
        -s $SAMPLENAME -e $ENVIRONMENT -q $QUEUEID -u $USER -p $PASSWORD
}

tmb_get_Tumor_Normal_Pair()
{
  getids_statement=" select samples.sampleID, samples.sampleName, sampleNormalPair.normalPairRunID, sampleNormalPair.normalSampleName from pipelineQueue \
  join samples on samples.sampleID = pipelineQueue.sampleID \
  join sampleNormalPair on sampleNormalPair.sampleID = pipelineQueue.sampleID where queueID='$QUEUEID';"

  while  read -r  sampleID sampleName normalPairRunID normalSampleName ; do

     TUMOR_SAMPLEID="$sampleID"
     TUMOR="$sampleName"
     NORMAL_RUNID="$normalPairRunID"
     NORMAL="$normalSampleName"

  done < <(mysql --host="$DB_HOST" --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$getids_statement" -N)

}

tmb_run_Alignment_paired()
{

  log_info "Running NORMAL/TUMOR Pair TMB Pipeline!"

  update_status "$QUEUEID" "Alignment" "$DB" "$USER"  "$PASSWORD" "$DB_HOST"

	/storage/apps/opt/parallel/bin/parallel --link "bash ${HOME_SHELLDIR}tmb_alignment_Interface.sh \
	    -s  {1} \
	    -f  {2} \
	    -o  $TMB_AL_OUT \
      -e  $ENVIRONMENT \
      -q  $QUEUEID  \
      -u  $USER \
      -p  $PASSWORD " echo ::: "$NORMAL" "$TUMOR" ::: "$FASTQ_DIR_NORMAL" "$FASTQ_DIR"

}


tmb_run_variant_caller_paired()
{
	if [ ! -f ${TMB_AL_OUT}${NORMAL}/Alignment/${NORMAL}.sort.rmdups.bam ] ; then
		log_info "Error: ${TMB_AL_OUT}${NORMAL}/Alignment/${NORMAL}.sort.rmdups.bam file not found"
		return
	fi

	if [ ! -f ${TMB_AL_OUT}${TUMOR}/Alignment/${TUMOR}.sort.rmdups.bam ] ; then
		log_info "Error: ${TMB_AL_OUT}${NORMAL}/Alignment/${NORMAL}.sort.rmdups.bam file not found"
		return
	fi

  update_status "$QUEUEID" "VariantCaller" "$DB" "$USER"  "$PASSWORD" "$DB_HOST"

	bash /storage/apps/pipelines/ngs_${ENVIRONMENT}/shell/tmb_VCaller_interface.sh \
	    -n ${TMB_AL_OUT}${NORMAL}/Alignment/${NORMAL}.sort.rmdups.bam \
	    -t ${TMB_AL_OUT}${TUMOR}/Alignment/${TUMOR}.sort.rmdups.bam \
	    -v varscan-strelka-mutect \
			-o ${TMB_VC_OUT} \
      -e $ENVIRONMENT \
      -i $TUMOR_SAMPLEID \
      -q $QUEUEID  \
      -u $USER \
      -p $PASSWORD

}

tmb_generate_stats()
{
  /storage/apps/opt/python3_4/bin/python3   ${HOME_PYTHONDIR}tmb_getSeqStat.py \
        "${TMB_AL_OUT}${TUMOR}/Alignment/"      "$TUMOR"  \
        "${TMB_AL_OUT}${NORMAL}/Alignment/"      "$NORMAL"  \
        "${TMB_VC_OUT}${TUMOR}_${NORMAL}/"


  wc -l /storage/analysis/environments/ngs_${ENVIRONMENT}/nextseqAnalysis/tmbAssay/${RUNNAME}/${TUMOR}/Single/${TUMOR}/Alignment/${TUMOR}.depth.filter2.exon_intersect.bed > ${TMB_VC_OUT}${TUMOR}_${NORMAL}/${TUMOR}.breadth_coverage

  /storage/apps/opt/python3_4/bin/python3  ${HOME_PYTHONDIR}tmb_vcQC.py \
  ${TMB_VC_OUT}${TUMOR}_${NORMAL}/${TUMOR}_${NORMAL}.variantcallers.combinev2.*.vep.parse.txt \
  ${TMB_VC_OUT}${TUMOR}_${NORMAL}/${TUMOR}_${NORMAL}.titv_ratio

  /storage/apps/opt/python3_4/bin/python3  ${HOME_PYTHONDIR}tmb_final_result.py \
  ${TMB_VC_OUT}${TUMOR}_${NORMAL}/${TUMOR}_${NORMAL}.seq_stats   \
  ${TMB_VC_OUT}${TUMOR}_${NORMAL}/${TUMOR}.breadth_coverage \
  ${TMB_VC_OUT}${TUMOR}_${NORMAL}/${TUMOR}_${NORMAL}_varscan_strelka_mutect_10_10_10.variants_results  \
  ${TMB_VC_OUT}${TUMOR}_${NORMAL}/${TUMOR}_${NORMAL}.tmb_result \
  ${TMB_VC_OUT}${TUMOR}_${NORMAL}/${TUMOR}_${NORMAL}.titv_ratio \
  ${TMB_VC_OUT}${TUMOR}_${NORMAL}/${TUMOR}_${NORMAL}.final_result.txt
}

tmb_run_dbUpdate()
{

    update_status "$QUEUEID" "UpdatingDatabase" "$DB" "$USER"  "$PASSWORD" "$DB_HOST"

    tmb_results="$TMB_VC_OUT${TUMOR}_${NORMAL}/${TUMOR}_${NORMAL}.tmb_result"

    tmb_results_statement="load data local infile '$tmb_results' into table sampleTumorMutationBurden FIELDS TERMINATED BY ',' (sampleID,TMBPair,TMBTotalVariants,TMBScore,TMBGroup)"
    mysql --host="$DB_HOST" --user="$USER" --password="$PASSWORD" --database="$DB" --execute="$tmb_results_statement"

    /storage/apps/opt/python3/bin/python3  ${HOME_PYTHONDIR}tmb_qc_plots.py \
    -i "$TUMOR_SAMPLEID" \
    -d "$DB_HOST" \
    -e "$ENVIRONMENT"  \
    -u "$USER"  \
    -p "$PASSWORD"

    update_status "$QUEUEID" "pipelineCompleted" "$DB" "$USER"  "$PASSWORD" "$DB_HOST"

}

# ##############################################################################
# main
# ##############################################################################
main()
{
    parse_options $*

    if [ $? -eq 0 ];then
			  echo "Error in previous step. Aborting $0"
        exit 0
    fi

		############################################################################
    # initialize variables
		############################################################################

    load_modules

    RUNFOLDER=$(ls -d /storage/instruments/$INSTRUMENT/*_"$RUNID"_*)
    RUNNAME=$(basename $RUNFOLDER)

    HOME_CODE="/storage/apps/pipelines/ngs_${ENVIRONMENT}/"
    HOME_RUNDIR="${HOME_CODE}run_files/"
    HOME_SHELLDIR="${HOME_CODE}shell/"
    HOME_PYTHONDIR="${HOME_CODE}python/"
    DB="ngs_${ENVIRONMENT}"
    DB_HOST="storage"

    HOME_ANALYSIS="/storage/analysis/environments/ngs_${ENVIRONMENT}/${INSTRUMENT}Analysis/${ASSAY}Assay/${RUNNAME}/${SAMPLENAME}/"


    update_status "$QUEUEID" "started" "$DB" "$USER"  "$PASSWORD" "$DB_HOST"


    if [ $ASSAY == "heme" ] ; then

      create_rundate

      exec >  >(tee -a ${HOME_ANALYSIS}process.log)
      exec 2> >(tee -a ${HOME_ANALYSIS}process.log >&2)

      show_pbsinfo


      VARIANT_VCF=""
      AMPLICON_FILE=""

  		create_dir ${HOME_ANALYSIS}variantCaller
  		create_dir ${HOME_ANALYSIS}variantAnalysis


      log_info " Running illumina Pipeline Interface BY qsub for :
  		assay : $ASSAY
  		instrument : $INSTRUMENT
  		runID : $RUNID
  		sampleName : $SAMPLENAME
  		environment : $ENVIRONMENT
  		queueID : $QUEUEID "

      heme_generate_variantFile

      heme_generate_ampliconFile

      log_info " amplicon file is - $AMPLICON_FILE"
      log_info " variant file is - $VARIANT_VCF"

      heme_pre_annotation

      heme_run_annotation

    elif [ $ASSAY == "tmb" ] ; then

      FASTQ_DIR="${RUNFOLDER}/out1/"
    	TMB_AL_OUT="${HOME_ANALYSIS}Single/"
    	TMB_VC_OUT="${HOME_ANALYSIS}Paired/"

      create_dir ${TMB_AL_OUT}
      create_dir ${TMB_VC_OUT}

      TUMOR_SAMPLEID=""
      TUMOR=""
      NORMAL_RUNID=""
      NORMAL=""

      tmb_get_Tumor_Normal_Pair

      FASTQ_DIR_NORMAL=""

      if [ "$NORMAL_RUNID" != "$RUNID" ] ; then

        NORMAL_RUNFOLDER=$(ls -d /storage/instruments/$INSTRUMENT/*_"$NORMAL_RUNID"_*)
        FASTQ_DIR_NORMAL="${NORMAL_RUNFOLDER}/out1/"

      else

        FASTQ_DIR_NORMAL="$FASTQ_DIR"

      fi

      LOG_FILE="${HOME_ANALYSIS}${TUMOR}_${NORMAL}.log"

      exec >  >(tee -a ${LOG_FILE})
      exec 2> >(tee -a ${LOG_FILE} >&2)

      show_pbsinfo

      log_info "Starting TMB Assay"
      ##############################################################################
      log_info "Parameters are-
      Tumor Sample-$TUMOR
      Normal Sample-$NORMAL
      Environment- $ENVIRONMENT
      Fastq Directory Tumor -$FASTQ_DIR
      Fastq Directory Normal-$FASTQ_DIR_NORMAL
      TMB Analysis Directory-$TMB_OUT
        $TMB_AL_OUT
        $TMB_VC_OUT
      Instrument-$INSTRUMENT
      Assay-$ASSAY
      code files:
      Home-$HOME_CODE
      Home Run Files-$HOME_RUNDIR
      Home Shell Files-$HOME_SHELLDIR
      Database-$DB"

    ##############################################################################

    # tmb_run_Alignment_paired

    # tmb_run_variant_caller_paired

    # tmb_generate_stats

    tmb_run_dbUpdate

    fi
}

# ##############################################################################
# run main
# ##############################################################################
main $*

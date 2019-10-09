#!/bin/bash

# ##############################################################################
# # functions
# ##############################################################################
display_usage()
{
cat <<EOF >> /dev/stderr

 USAGE: $0

 OPTIONS:
 	-s SAMPLE NAME
 	-f FASTQ Directory
	-o Output Directory

EOF
}

parse_options()
{
    IMPORT=0

		while getopts "hz:s:f:o:e:q:u:p:" opt; do

				case $opt in
					h)
					display_usage
					exit 1
					;;
	        z)
					IMPORT=$OPTARG
					;;
					s)
					SAMPLE=$OPTARG
					;;
					f)
					FASTQ_DIR=$OPTARG
					;;
					o)
					OUTPUT_DIR=$OPTARG
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
  source /home/pipelines/ngs_${ENVIRONMENT}/shell/modules/ngs_align.sh
  source /home/pipelines/ngs_${ENVIRONMENT}/shell/modules/ngs_varscan.sh
  source /home/pipelines/ngs_${ENVIRONMENT}/shell/modules/ngs_vep.sh

}

run_trimmomatic()
{
	log_info "Running Trimmomatic: Removing sequences < Q20 sample- $SAMPLE"
	trimmomatic="java -jar /opt/trimmomatic/Trimmomatic-0.33/trimmomatic-0.33.jar PE -phred33 -threads 8 \
	${FASTQ_DIR}${SAMPLE}_S*_R1_001.fastq.gz \
	${FASTQ_DIR}${SAMPLE}_S*_R2_001.fastq.gz \
	${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}_filt_paired_R1_001.fastq.gz \
	${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}_filt_unpaired_R1_001.fastq.gz \
	${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}_filt_paired_R2_001.fastq.gz \
	${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}_filt_unpaired_R2_001.fastq.gz \
	TRAILING:20 \
	AVGQUAL:20 \
	SLIDINGWINDOW:10:20 \
	MINLEN:30 \
  HEADCROP:20"
	($trimmomatic) 2>&1 | tee ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.trimmomatic.summary.txt
}

run_alignment()
{
	log_info "Running bwa mem aligner: $SAMPLE"

  tmb_bwaAlign $SAMPLE  \
	$REF_GENOME  \
	$MAP_QUALITY \
	${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}_filt_paired_R1_001.fastq.gz \
	${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}_filt_paired_R2_001.fastq.gz  \
	$OUTPUT_DIR_SAMPLE_ALIGNMENT  \
	$LOG_FILE
}

sort_bam()
{
log_info "Generating sorted bam by coordinate sample- $SAMPLE"
java -jar /opt/picard2/picard.jar SortSam \
          I=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.bam  \
          O=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.sorted.bam  \
          SORT_ORDER=coordinate
}

generate_alignStat()
{
log_info "Generating alignment stat sample- $SAMPLE"
java -jar /opt/picard2/picard.jar CollectAlignmentSummaryMetrics \
          R=$REF_GENOME \
          I=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.sorted.bam \
          O=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.sorted.bam.alignmentMetrics.txt
}

# remove_duplicates()
# {
# log_info "Removing duplicates sample- $SAMPLE "
# java -jar /opt/picard2/picard.jar MarkDuplicates \
#           REMOVE_DUPLICATES=true \
#           INPUT=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.sorted.bam \
#           OUTPUT=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.sorted.rmdups.bam \
#           METRICS_FILE=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.sorted.rmdups.bam.metrics.txt
# }

generate_bamIndex()
{
log_info "Generating bam index sample- $SAMPLE"
java -jar /opt/picard2/picard.jar BuildBamIndex \
          I=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.sorted.bam
}


calculate_targetCoverage()
{
log_info "Generating CalculateHsMetrics sample- $SAMPLE "
# # java -jar /opt/picard/picard-tools-1.134/picard.jar BedToIntervalList  I=myeloid_design.bed O=myeloid_design.interval_list SD=/doc/ref/ref_genome/ucsc.hg19.dict
java -jar /opt/picard/picard-tools-1.134/picard.jar CalculateHsMetrics \
          I=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.sorted.bam  \
          O=${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.output_hs_metrics.txt \
          R=$REF_GENOME \
          BAIT_INTERVALS= /home/environments/ngs_test/assayCommonFiles/hemeAssay/myeloid_design.interval_list \
          TARGET_INTERVALS= /home/environments/ngs_test/assayCommonFiles/hemeAssay/myeloid_design.interval_list
}

heme_generate_variantFile()
{

  heme_varscan ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.sorted.bam  ${OUTPUT_DIR_SAMPLE_ALIGNMENT}


  bash ${HOME_SHELL}heme_parseVarScan.sh \
          -s ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.snp.txt \
          -i ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.indel.txt \
          -o ${OUTPUT_DIR_SAMPLE_ALIGNMENT}  \
          -e $ENVIRONMENT

   VARIANT_VCF=$(ls ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.comb.vcf)
}


heme_pre_annotation()
{

  HEME_EXCLUDED_DESIGN="/home/doc/ref/Heme/trusight-myeloid-amplicon-track.excluded.bed"
  #bedtools -u flag to write original A entry once if any overlaps found in B
  /opt/software/bedtools-2.17.0/bin/bedtools intersect -u -a $VARIANT_VCF -b $HEME_EXCLUDED_DESIGN > ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.filter.vcf


  # HEME_EXCLUDED_AMPLICON="/home/doc/ref/Heme/excludedAmplicons.txt"
  # # -v means "invert the match" in grep, in other words, return all non matching lines.
  # grep -v -f $HEME_EXCLUDED_AMPLICON $AMPLICON_FILE > ${HOME_ANALYSIS}variantAnalysis/${SAMPLENAME}.amplicon.filter.txt

}

run_vep()
{
	log_info "Running VEP"

  # vep_83 "${HOME}variantAnalysis/${SAMPLENAME}.filter.vcf"  "${HOME}variantAnalysis/${SAMPLENAME}.filter.vep.vcf"

  vep_94_panel "${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.filter.vcf"  "${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.filter.vep.vcf"

  log_info "Completed VEP"
}

parse_vep()
{
  log_info "Parse VEP"

	python ${HOME_PYTHON}parseVEP_v2.py \
					parseIlluminaNextseq   \
					-I ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.filter.vep.vcf\
					-o ${OUTPUT_DIR_SAMPLE_ALIGNMENT}${SAMPLE}.filter.vep.parse.vcf
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

	DIR_SCRIPT="/home/pipelines/ngs_${ENVIRONMENT}/"
  HOME_PYTHON="/home/pipelines/ngs_${ENVIRONMENT}/python/"
  HOME_SHELL="/home/pipelines/ngs_${ENVIRONMENT}/shell/"
	REF_GENOME="/home/doc/ref/ref_genome/ucsc.hg19.fasta"
	MAP_QUALITY="20"

  load_modules

	OUTPUT_DIR_SAMPLE="${OUTPUT_DIR}${SAMPLE}/"
	create_dir $OUTPUT_DIR_SAMPLE

	LOG_FILE="${OUTPUT_DIR}${SAMPLE}.log"


  exec >  >(tee -a ${LOG_FILE})
  exec 2> >(tee -a ${LOG_FILE} >&2)

  log_info "Starting alignment"

	OUTPUT_DIR_SAMPLE_ALIGNMENT="${OUTPUT_DIR_SAMPLE}Alignment/"
	create_dir $OUTPUT_DIR_SAMPLE_ALIGNMENT


	log_info "
	Starting - heme_Interface.sh
	SAMPLE - $SAMPLE
	FASTQ_DIR - $FASTQ_DIR
	OUTPUT_DIR - $OUTPUT_DIR
	OUTPUT_DIR_SAMPLE - $OUTPUT_DIR_SAMPLE
	OUTPUT_DIR_SAMPLE_ALIGNMENT - $OUTPUT_DIR_SAMPLE_ALIGNMENT
	LOG_FILE - $LOG_FILE "

  # run_trimmomatic
  #
  # run_alignment
  #
  # sort_bam
  #
  # generate_alignStat
  #
  # #remove_duplicates
  #
  # generate_bamIndex
  #
  # calculate_targetCoverage


  VARIANT_VCF=""


  heme_generate_variantFile

  heme_pre_annotation

  run_vep

  parse_vep


}

# ##############################################################################
# run main
# ##############################################################################
main $*

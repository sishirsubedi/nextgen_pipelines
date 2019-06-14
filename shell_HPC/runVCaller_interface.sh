#!/bin/bash

# ##############################################################################
# # functions
# ##############################################################################

display_usage()
{
cat <<EOF >> /dev/stderr

 USAGE: $0

 OPTIONS:
 	-n Normal Sample
 	-t Tumor Sample
	-v Variant Callers
	-o Output Directory

EOF
}

parse_options()
{
    IMPORT=0

		while getopts "hz:n:t:v:o:" opt; do

				case $opt in
					h)
					display_usage
					exit 1
					;;
	        z)
					IMPORT=$OPTARG
					;;
					n)
					NORMAL_BAM=$OPTARG
					;;
					t)
					TUMOR_BAM=$OPTARG
					;;
					v)
					VARIANT_CALLERS=$OPTARG
					;;
					o)
					OUT_DIR=$OPTARG
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
	source /home/pipelines/ngs_${ENVIRONMENT}/shell/modules/ngs_varscan.sh
}

run_UNpairedVC()
{
for v_caller in $VCS;do

	create_dir ${SAMPLE_DIR}/${v_caller}

	if [[ "$v_caller" = "varscan" ]];then

		echo "Starting : " $v_caller
		heme_varscan $NORMAL_BAM   ${SAMPLE_DIR}/${v_caller}/

	elif [[ "$v_caller" = "strelka" ]]; then

		echo "Starting : " $v_caller

	elif [[ "$v_caller" = "mutect" ]]; then

		echo "Starting : " $v_caller

	fi

done
}

run_pairedVC()
{

	DEPTH="10"
	NALF="10"
	TALF="2"

	for v_caller in $VCS;do

		create_dir ${SAMPLE_DIR}/${v_caller}

		if [[ "$v_caller" = "varscan" ]];then
			echo "Starting : " $v_caller
			bash ${DIR_SCRIPT}shell/runVCaller_varscan.sh $SAMPLE $REF_GENOME_1  $NORMAL_BAM  $TUMOR_BAM  ${SAMPLE_DIR}/${v_caller}/  $ENVIRONMENT  $DEPTH  $NALF  $TALF
		elif [[ "$v_caller" = "strelka" ]]; then
			echo "Starting : " $v_caller
			bash ${DIR_SCRIPT}shell/runVCaller_strelka.sh $SAMPLE $REF_GENOME_1  $NORMAL_BAM  $TUMOR_BAM  ${SAMPLE_DIR}/${v_caller}/  $ENVIRONMENT $DEPTH $NALF $TALF
		elif [[ "$v_caller" = "mutect" ]]; then
			echo "Starting : " $v_caller
			bash ${DIR_SCRIPT}shell/runVCaller_mutect.sh $SAMPLE $REF_GENOME_1  $NORMAL_BAM  $TUMOR_BAM  ${SAMPLE_DIR}/${v_caller}/  $ENVIRONMENT $DEPTH $NALF $TALF
	  fi

	done


	/opt/python3/bin/python3 /home/hhadmin/scripts/bioinfoTools/03_variant_report/3_compare_samples_3Venn.py   \
	${SAMPLE_DIR}/varscan/${SAMPLE}.varscan.${DEPTH}_${NALF}_${TALF}  \
	${SAMPLE_DIR}/strelka/${SAMPLE}.strelka.${DEPTH}_${NALF}_${TALF} \
	${SAMPLE_DIR}/mutect/${SAMPLE}.mutect.${DEPTH}_${NALF}_${TALF} \
	${SAMPLE_DIR}/ "${DEPTH}_${NALF}_${TALF}"


	##### combine vcfs from all variant callers
	/opt/python3/bin/python3 ${DIR_SCRIPT}python/combineVCFs.py  "$SAMPLE"  "$SAMPLE_DIR"  "$ENVIRONMENT"  "${DEPTH}_${NALF}_${TALF}"


	tail -n +2 "${SAMPLE_DIR}/${SAMPLE}.variantcallers.combine.${DEPTH}_${NALF}_${TALF}" > ${SAMPLE_DIR}/${SAMPLE}.variantcallers.combinev2.${DEPTH}_${NALF}_${TALF}
	echo " $currentdate    INFO  -  running VEP"
	/opt/vep_94/ensembl-tools-release-94/vep_94/ensembl-vep/vep \
	-i ${SAMPLE_DIR}/${SAMPLE}.variantcallers.combinev2.${DEPTH}_${NALF}_${TALF} \
	-o ${SAMPLE_DIR}/${SAMPLE}.variantcallers.combinev2.${DEPTH}_${NALF}_${TALF}.vep \
	--offline \
	--dir_cache /opt/vep_94/ensembl-tools-release-94/cache \
	--vcf \
	--refseq \
	--pick_allele \
	--sift p \
	--polyphen p \
	--hgvs \
	--symbol \
	--vcf \
	--pubmed \
	--fasta $REF_GENOME_2 \
	--force_overwrite


	# ##### combine vcfs from all variant callers
	/opt/python3/bin/python3 ${DIR_SCRIPT}python/parseVEP_exome.py  "$SAMPLE"  "$SAMPLE_DIR"  "$ENVIRONMENT"  "${DEPTH}_${NALF}_${TALF}"

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


	ENVIRONMENT="test"
	DIR_SCRIPT="/home/pipelines/ngs_${ENVIRONMENT}/"
	REF_GENOME_1="/home/doc/ref/ref_genome/ucsc.hg19.fasta"
	REF_GENOME_2="/home/doc/ref/ref_genome/Homo_sapiens.GRCh37.73.dna_sm.primary_assembly.fa"

  load_modules
  log_info "Starting"

	if [[ $TUMOR_BAM = "NONE" ]];then
		SAMPLE=$(echo "${NORMAL_BAM##*/}"| tr "." "\n" | head -1 )
	else
		SAMPLE=$(echo "${TUMOR_BAM##*/}"| tr "." "\n" | head -1 )
		SAMPLE+="_"$(echo "${NORMAL_BAM##*/}"| tr "." "\n" | head -1 )
	fi


	SAMPLE_DIR=${OUT_DIR}${SAMPLE}
	create_dir $SAMPLE_DIR

	#get all variant callers in list
	VCS=$(echo $VARIANT_CALLERS | tr "-" "\n")

	if [[ $TUMOR_BAM = "NONE" ]];then
		run_UNpairedVC
	else
		run_pairedVC
	fi

}

# ##############################################################################
# run main
# ##############################################################################
main $*

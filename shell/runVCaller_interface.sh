#!/usr/bin/env bash

ENV="test"
REF_1="/home/doc/ref/ref_genome/ucsc.hg19.fasta"
#################################################
# Parsing arguments
#################################################

if [ "$#" -eq 0 ]; then
echo "Usage: runVariantCaller.sh"
echo "-n Normal BAM"
echo "-t TUMOR BAM"
echo "-v List of Variant Callers"
echo "-o Output Dir"
exit
fi

while getopts :n:t:v:o: option; do
	case "$option" in
    n) NORMAL_BAM="$OPTARG" ;;
    t) TUMOR_BAM="$OPTARG" ;;
    v) VARIANT_CALLERS="$OPTARG" ;;
    o) OUT_DIR="$OPTARG" ;;
    :) echo "Option -$OPTARG requires an argument." ;;
	  \?) echo "Invalid option: -$OPTARG" ;;
	esac
done

#get sample name
SAMPLE=$(echo "${TUMOR_BAM##*/}"| tr "." "\n" | head -1 )
SAMPLE+="_"$(echo "${NORMAL_BAM##*/}"| tr "." "\n" | head -1 )

SAMPLE_DIR=${OUT_DIR}${SAMPLE}

if [ ! -d $SAMPLE_DIR ] ; then
	sudo mkdir $SAMPLE_DIR
	sudo chmod 777 $SAMPLE_DIR
fi

# get all variant callers in list
VCS=$(echo $VARIANT_CALLERS | tr "-" "\n")

for v_caller in $VCS;do

	if [ ! -d ${SAMPLE_DIR}/${v_caller} ] ; then
		sudo mkdir ${SAMPLE_DIR}/${v_caller}
		sudo chmod 777 ${SAMPLE_DIR}/${v_caller}
	fi


	if [[ "$v_caller" = "varscan" ]];then
		echo "Starting : " $v_caller
		# bash /var/pipelines_ngs_${ENV}/shell/runVCaller_varscan.sh $SAMPLE $REF_1  $NORMAL_BAM  $TUMOR_BAM  ${SAMPLE_DIR}/${v_caller}/  $ENV
	elif [[ "$v_caller" = "mutect" ]]; then
		echo "Starting : " $v_caller
		# bash /var/pipelines_ngs_${ENV}/shell/runVCaller_mutect.sh $REF_1  $NORMAL_BAM  $TUMOR_BAM  ${SAMPLE_DIR}/${v_caller}/  $ENV
	elif [[ "$v_caller" = "strelka" ]]; then
		echo "Starting : " $v_caller
		/var/pipelines_ngs_${ENV}/shell/runVCaller_strelka.sh $REF_1  $NORMAL_BAM  $TUMOR_BAM  ${SAMPLE_DIR}/${v_caller}/  $ENV
	elif [[ "$v_caller" = "ebcall" ]]; then
		echo "Starting : " $v_caller
		# /var/pipelines_ngs_${ENV}/shell/runVCaller_ebcall.sh $REF_1  $NORMAL_BAM  $TUMOR_BAM  ${SAMPLE_DIR}/${v_caller}/   $ENV
  fi
done




#
#
# echo " $currentdate    INFO  -  running VEP"
# /opt/vep_94/ensembl-tools-release-94/vep_94/ensembl-vep/vep \
# -i /home/hhadmin/exome_pipeline/02_variantCalling/$sample/"$sample"_snp_indel_filter.vcf \
# -o /home/hhadmin/exome_pipeline/02_variantCalling/$sample/"$sample"_snp_indel_filter.vep.vcf.txt \
# --offline \
# --dir_cache /opt/vep_94/ensembl-tools-release-94/cache \
# --vcf \
# --refseq \
# --pick_allele \
# --sift p \
# --polyphen p \
# --hgvs \
# --symbol \
# --vcf \
# --pubmed \
# --fasta /home/doc/ref/ref_genome/Homo_sapiens.GRCh37.73.dna_sm.primary_assembly.fa \
# --force_overwrite

#!/usr/bin/env bash

ENV="test"
DIR_SCRIPT="/home/pipelines/ngs_${ENV}/"
REF_GENOME_1="/home/doc/ref/ref_genome/ucsc.hg19.fasta"
REF_GENOME_2="/home/doc/ref/ref_genome/Homo_sapiens.GRCh37.73.dna_sm.primary_assembly.fa"


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
	mkdir $SAMPLE_DIR
	chmod 775 $SAMPLE_DIR
fi

#get all variant callers in list
VCS=$(echo $VARIANT_CALLERS | tr "-" "\n")

DEPTH="0"
NALF="100"
TALF="0"

for v_caller in $VCS;do

	if [ ! -d ${SAMPLE_DIR}/${v_caller} ] ; then
		mkdir ${SAMPLE_DIR}/${v_caller}
		chmod 775 ${SAMPLE_DIR}/${v_caller}
	fi

	if [[ "$v_caller" = "varscan" ]];then
		echo "Starting : " $v_caller
		bash ${DIR_SCRIPT}shell/runVCaller_varscan.sh $SAMPLE $REF_GENOME_1  $NORMAL_BAM  $TUMOR_BAM  ${SAMPLE_DIR}/${v_caller}/  $ENV  $DEPTH  $NALF  $TALF
	elif [[ "$v_caller" = "strelka" ]]; then
		echo "Starting : " $v_caller
		bash ${DIR_SCRIPT}shell/runVCaller_strelka.sh $SAMPLE $REF_GENOME_1  $NORMAL_BAM  $TUMOR_BAM  ${SAMPLE_DIR}/${v_caller}/  $ENV $DEPTH $NALF $TALF
	elif [[ "$v_caller" = "mutect" ]]; then
		echo "Starting : " $v_caller
		bash ${DIR_SCRIPT}shell/runVCaller_mutect.sh $SAMPLE $REF_GENOME_1  $NORMAL_BAM  $TUMOR_BAM  ${SAMPLE_DIR}/${v_caller}/  $ENV $DEPTH $NALF $TALF
  fi

done



/opt/python3/bin/python3 /home/hhadmin/scripts/bioinfoTools/03_variant_report/3_compare_samples_3Venn.py   \
${SAMPLE_DIR}/varscan/${SAMPLE}.varscan.${DEPTH}_${NALF}_${TALF}  \
${SAMPLE_DIR}/strelka/${SAMPLE}.strelka.${DEPTH}_${NALF}_${TALF} \
${SAMPLE_DIR}/mutect/${SAMPLE}.mutect.${DEPTH}_${NALF}_${TALF} \
${SAMPLE_DIR}/ "${DEPTH}_${NALF}_${TALF}"


##### combine vcfs from all variant callers
/opt/python3/bin/python3 ${DIR_SCRIPT}python/combineVCFs.py  "$SAMPLE"  "$SAMPLE_DIR"  "$ENV"  "${DEPTH}_${NALF}_${TALF}"


# tail -n +2 "${SAMPLE_DIR}/${SAMPLE}.variantcallers.combine.${DEPTH}_${NALF}_${TALF}" > ${SAMPLE_DIR}/${SAMPLE}.variantcallers.combinev2.${DEPTH}_${NALF}_${TALF}
# echo " $currentdate    INFO  -  running VEP"
# /opt/vep_94/ensembl-tools-release-94/vep_94/ensembl-vep/vep \
# -i ${SAMPLE_DIR}/${SAMPLE}.variantcallers.combinev2.${DEPTH}_${NALF}_${TALF} \
# -o ${SAMPLE_DIR}/${SAMPLE}.variantcallers.combinev2.${DEPTH}_${NALF}_${TALF}.vep \
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
# --fasta $REF_GENOME_2 \
# --force_overwrite
#
#
#
#
# # ##### combine vcfs from all variant callers
# /opt/python3/bin/python3 ${DIR_SCRIPT}python/parseVEP_exome.py  "$SAMPLE"  "$SAMPLE_DIR"  "$ENV"  "${DEPTH}_${NALF}_${TALF}"

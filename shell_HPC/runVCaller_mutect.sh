#!/usr/bin/env bash

#################################################
# Mutect
#################################################

SAMPLE=$1
REF=$2
NORMAL_BAM=$3
TUMOR_BAM=$4
OUT_DIR=$5
ENV=$6

echo " starting mutect:
      REF : $REF
      NORMAL_BAM : $NORMAL_BAM
      TUMOR_BAM : $TUMOR_BAM
      OUT_DIR : $OUT_DIR "


suffix=".bam"
NORMAL_BAM_CHR=$(echo $NORMAL_BAM | sed -e "s/$suffix$//" )
TUMOR_BAM_CHR=$(echo $TUMOR_BAM | sed -e "s/$suffix$//" )

# create chromosome bam using bamtools
# normal_split="/opt/bamtools/bamtools/bin/bamtools split -in  $NORMAL_BAM -reference"
# tumor_split="/opt/bamtools/bamtools/bin/bamtools split -in  $TUMOR_BAM -reference"
# (echo $tumor_split; echo $normal_split) | /opt/parallel/bin/parallel


# index all chr specific bam files
# for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
#
#   echo ${NORMAL_BAM_CHR}.REF_chr${chr}.bam
#   echo ${TUMOR_BAM_CHR}.REF_chr${chr}.bam
#
#   normal_index="/opt/samtools19/bin/samtools index ${NORMAL_BAM_CHR}.REF_chr${chr}.bam"
#   tumor_index="/opt/samtools19/bin/samtools index ${TUMOR_BAM_CHR}.REF_chr${chr}.bam"
#
#   (echo $tumor_index; echo $normal_index) | /opt/parallel/bin/parallel
#
# done

#
# for counter in 1 6 11 16 21; do
#   ##check chromosome input file
#     if [ -f  ${OUT_DIR}active_chromosomes.txt ]; then
#        rm -f ${OUT_DIR}active_chromosomes.txt
#     fi
#
#   ## write to input file
#     if [ $counter = "21" ]; then
#       for inner_counter in 21 22 X Y;do
#         chr=$inner_counter
#         echo $chr >> ${OUT_DIR}active_chromosomes.txt
#       done
#     else
#       for inner_counter in 0 1 2 3 4;do
#         chr=$((counter + inner_counter))
#         echo $chr >> ${OUT_DIR}active_chromosomes.txt
#       done
#     fi
#
#   cat ${OUT_DIR}active_chromosomes.txt
#
#
#   ## run gatk in parallel
#   cat ${OUT_DIR}active_chromosomes.txt | /opt/parallel/bin/parallel "java -jar /opt/GATK4/GenomeAnalysisTK.jar \
#                                       -T MuTect2 \
#                                       -R $REF \
#                                       -I:tumor ${TUMOR_BAM_CHR}.REF_chr{}.bam \
#                                       -I:normal ${NORMAL_BAM_CHR}.REF_chr{}.bam \
#                                       -o ${OUT_DIR}chr{}_gatk.output.vcf \
#                                       -L chr{} \
#                                       --min_base_quality_score 30 "
# done
#



for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
  java -jar /opt/GATK4/GenomeAnalysisTK.jar \
       -R $REF \
       -T VariantsToTable \
       -V ${OUT_DIR}chr"$chr"_gatk.output.vcf \
       -F CHROM -F POS -F REF -F ALT -F FILTER \
       -F ECNT -F HCNT -F MAX_ED -F MIN_ED -F NLOD -F RPA -F RU=CA -F STR -F TLOD \
       -GF GT -GF AD -GF AF -GF ALT_F1R2 -GF ALT_F2R1 -GF FOXOG -GF PGT  -GF PID  -GF QSS  -GF REF_F1R2 -GF REF_F2R1 \
       -o ${OUT_DIR}chr"$chr"_gatk.output.vcf.filter.txt
done
#
/opt/python3/bin/python3 /home/hhadmin/pipelines_ngs_${ENV}/python/parseGATK.py  "$SAMPLE"  "$OUT_DIR"  "$ENV"
# /opt/python3/bin/python3 /var/pipelines_ngs_${ENV}/python/parseGATK_paraEST.py  "$SAMPLE"  "$OUT_DIR"  "$ENV"

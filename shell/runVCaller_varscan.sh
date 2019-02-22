#!/usr/bin/env bash

#################################################
# Varscan
#################################################
SAMPLE=$1
REF=$2
NORMAL_BAM=$3
TUMOR_BAM=$4
OUT_DIR=$5
ENV=$6


echo " starting varscan:
REF : $REF
NORMAL_BAM :$NORMAL_BAM
TUMOR_BAM :$TUMOR_BAM
OUT_DIR :$OUT_DIR
"

#
# tumor_pileup="/opt/samtools19/bin/samtools  mpileup  -f $REF  $TUMOR_BAM > ${OUT_DIR}${SAMPLE}.tumor.pileup"
# normal_pileup="/opt/samtools19/bin/samtools  mpileup -f $REF  $NORMAL_BAM > ${OUT_DIR}${SAMPLE}.normal.pileup"
#
# (echo $tumor_pileup; echo $normal_pileup) | /opt/parallel/bin/parallel
#
# java -jar /opt/varscan/VarScan.v2.3.9.jar somatic \
#           ${OUT_DIR}${SAMPLE}.normal.pileup \
#           ${OUT_DIR}${SAMPLE}.tumor.pileup \
#           ${OUT_DIR}${SAMPLE}.varscan.output \
#           --min-avg-qual 30

/opt/python3/bin/python3 /home/hhadmin/pipelines_ngs_${ENV}/python/parseVarscan.py  $SAMPLE  $OUT_DIR  $ENV
# /opt/python3/bin/python3 /var/pipelines_ngs_${ENV}/python/parseVarscan_paraEst.py  $SAMPLE  $OUT_DIR  $ENV

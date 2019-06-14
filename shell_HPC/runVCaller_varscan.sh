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
DEPTH=$7
NALF=$8
TALF=$9


echo " starting varscan:
REF : $REF
NORMAL_BAM :$NORMAL_BAM
TUMOR_BAM :$TUMOR_BAM
OUT_DIR :$OUT_DIR
"

java -jar /opt/varscan/VarScan.v2.3.9.jar somatic \
          <(/opt/samtools19/bin/samtools  mpileup  -f $REF  $NORMAL_BAM ) \
          <(/opt/samtools19/bin/samtools  mpileup  -f $REF  $TUMOR_BAM) \
          ${OUT_DIR}${SAMPLE}.varscan.output \
          --min-avg-qual 30


/opt/python3/bin/python3 /home/pipelines/ngs_${ENV}/python/parseVarscan.py  "$SAMPLE"  "$OUT_DIR"  "$ENV" "$DEPTH" "$NALF" "$TALF"



# /opt/python3/bin/python3 /home/pipelines/ngs_${ENV}/python/filterCreDesign.py  ${OUT_DIR}  ${SAMPLE}.varscan.${DEPTH}_${NALF}_${TALF}  $ENV

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

/storage/apps/opt/java/jdk1.8.0_191/bin/java -jar /storage/apps/opt/varscan/VarScan.v2.3.9.jar somatic \
          <(/storage/apps/opt/samtools/bin/samtools  mpileup  -f $REF  $NORMAL_BAM ) \
          <(/storage/apps/opt/samtools/bin/samtools  mpileup  -f $REF  $TUMOR_BAM) \
          ${OUT_DIR}${SAMPLE}.varscan.output \
          --min-avg-qual 30 \
          --min-var-freq 0.1

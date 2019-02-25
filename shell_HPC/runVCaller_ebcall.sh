#!/usr/bin/env bash

#################################################
# Varscan
#################################################

REF=$1
NORMAL_BAM=$2
TUMOR_BAM=$3
OUT_DIR=$4

echo " starting ebcall:
REF : $REF
NORMAL_BAM : $NORMAL_BAM
TUMOR_BAM : $TUMOR_BAM
OUT_DIR : $OUT_DIR
"

sh ebCall_v2.sh $TUMOR_BAM $NORMAL_BAM $OUT_DIR /home/hhadmin/exome_pipeline/02_variantCalling/COV_1_COV_2/ebcall/gen_ref.txt

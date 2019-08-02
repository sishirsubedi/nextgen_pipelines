#!/bin/bash
input="/home/environments/ngs_test/nextseqAnalysis/exomeAnalysis/varscan_samples.csv"
OLDIFS=$IFS
IFS=","
while read run sample
do
  sampleName="$(echo -e "${sample}" | tr -d '[:space:]')"
  # echo "--$run--$sampleName--"
  qsub -F "$run  $sampleName" /home/pipelines/ngs_test/shell/tmb_var_Interface.sh
  # bash /home/pipelines/ngs_test/shell/tmb_var_Interface.sh  $run  $sampleName

done < "$input"
IFS=$OLDIFS

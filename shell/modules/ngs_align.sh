#!/bin/bash


heme_bwaAlign()
{
  fastq1=$1
  fastq2=$2
  outDir=$3
  fileName=${fastq1##*/}
  sample=${fileName%%_*}


  echo "#######Aligning " $sample
  /opt/bwa/bwa-0.7.12/bwa mem -t 8 /home/doc/ref/ref_genome/ucsc.hg19.fasta "$fastq1" "$fastq2" |/opt/samtools/samtools view -bS - > "$outDir"/"$sample".bam


  echo "#####generate stat " $sample
  /opt/samtools/samtools flagstat "$outDir"/"$sample".bam > "$outDir"/"$sample".flagstat


  echo "####sorting bam " $sample
  ##sort##
  java -Xmx1g -jar /opt/picard/picard-tools-1.134/picard.jar AddOrReplaceReadGroups \
  INPUT="$outDir"/"$sample".bam \
  OUTPUT="$outDir"/"$sample".sort.bam \
  SORT_ORDER=coordinate \
  RGID="$sample" \
  RGLB=1 \
  RGPL=illumina \
  RGPU=1 \
  RGSM="$sample"


  ##index###
  echo "#####Indexing " $sample
  /storage/apps/opt/samtools/samtools index "$outDir"/"$sample".sort.bam
}


tmb_bwaAlign()
{

  SAMPLE=$1
  REF=$2
  MAP_QUALITY=$3
  FASTQ1=$4
  FASTQ2=$5
  OUT_DIR=$6
  LOG_FILE=$7


  echo "Starting BWA alignment and generating bam file for:
  SAMPLE - $SAMPLE
  REF - $REF
  FASTQ1 - $FASTQ1
  FASTQ2 - $FASTQ2
  OUT_DIR - $OUT_DIR
  LOG_FILE - $LOG_FILE
  "
  /storage/apps/opt/bwa/bwa-0.7.12/bwa mem -M -t 32 -R "@RG\tID:${SAMPLE}\tLB:1\tSM:${SAMPLE}\tPL:ILLUMINA\tPU:1"  "$REF" "$FASTQ1" "$FASTQ2" | /storage/apps/opt/samtools/bin/samtools view --threads 32 -bf 0x2 -q "$MAP_QUALITY" > ${OUT_DIR}${SAMPLE}.bam


}

SAMPLE=$1
REF=$2
FASTQ1=$3
FASTQ2=$4
OUT_DIR=$5
LOG_FILE=$6


function log {
 MESSAGE=$1
 TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
 SCRIPT="runAlignment_Interface.sh"
 echo " [ $TIMESTAMP ] [ $SCRIPT ] : $MESSAGE "
 echo " [ $TIMESTAMP ] [ $SCRIPT ] : $MESSAGE " >>${LOG_FILE}
}

log "Starting BWA alignment and generating bam file for:
SAMPLE - $SAMPLE
REF - $REF
FASTQ1 - $FASTQ1
FASTQ2 - $FASTQ2
OUT_DIR - $OUT_DIR
LOG_FILE - $LOG_FILE
"
/opt/bwa/bwa-0.7.12/bwa mem -M -t 16 -R "@RG\tID:${SAMPLE}\tLB:1\tSM:${SAMPLE}\tPL:ILLUMINA\tPU:1"  "$REF" "$FASTQ1" "$FASTQ2" | /opt/samtools19/bin/samtools view --threads 16 -bf 0x2 -q 20 > ${OUT_DIR}${SAMPLE}.bam

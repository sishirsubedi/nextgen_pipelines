fastq1=$1
fastq2=$2
outDir=$3
fileName=${fastq1##*/}
sample=${fileName%%_*}
echo "#######Aligning " $sample
/opt/bwa/bwa-0.7.12/bwa mem /home/doc/ref/ref_genome/ucsc.hg19.fasta "$fastq1" "$fastq2" |/opt/samtools/samtools view -bS - > "$outDir"/"$sample".bam
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
/opt/samtools/samtools index "$outDir"/"$sample".sort.bam
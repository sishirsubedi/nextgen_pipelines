bam=$1
outDir=$2
fileName=${bam##*/}
sample=${fileName%%.*}

##Varscan call SNPs##
/opt/samtools/samtools mpileup -f /home/doc/ref/ref_genome/ucsc.hg19.fasta "$bam" |\
java -jar /opt/varscan/VarScan.v2.3.9.jar pileup2snp > "$outDir"/"$sample".snp.txt

##Varscan call indels##
/opt/samtools/samtools mpileup -f /home/doc/ref/ref_genome/ucsc.hg19.fasta "$bam" |\
java -jar /opt/varscan/VarScan.v2.3.9.jar pileup2indel > "$outDir"/"$sample".indel.txt


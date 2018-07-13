#!/bin/bash

#hotspot=$1
#outputDir=$2
#outPrefix=$3
#bam=$4


if [ $# -eq 0 ]
then
echo "Usage: readCount.sh"
echo "-l hotspot file, first 4 columns: chr, pos, ref, alt"
echo "-d output directory"
echo "-o output prefix" 
echo "-b bamfile"
exit
fi


if test $# -gt 0
	then
	while getopts :l:d:o:b: opt
	do
	case $opt in
	
	l)
		hotspot=$OPTARG
		;;
	d)
		outputDir=$OPTARG
		;;
	o)
		outPrefix=$OPTARG
		;;
	b)
		bam=$OPTARG
		;;
	:)
		echo "Option -$OPTARG requires an argument."
		;;
	\?)
		echo "Invalid option: -$OPTARG"
	esac
	done
	shift $((OPTIND-1))
fi


/opt/samtools-1.3/samtools-1.3/bin/samtools mpileup \
-B \
-C 0 \
-t INFO/AD,INFO/ADF,INFO/ADR \
-f /doc/ref/ref_genome/ucsc.hg19.fasta \
-l $hotspot \
-o $outputDir/$outPrefix.pileup.vcf  \
-v \
-A \
-Q 0 \
-u \
-x \
-d 100000 \
-F 0 \
-L 100000 \
$bam

python ~/code/python/parsePileup.py -I $outputDir/$outPrefix.pileup.vcf -o $outputDir/$outPrefix.pileup.txt

join -t$'\t' -a 1 -1 1 -2 1 \
<(awk -v OFS='\t' -F'\t' '{print $1"|"$2, $3,$4, $5,$6,$7,$8,$9,$10}' /projects/simulation/readCount/lungHotspot.linux.txt |sort -t$'\t') \
<(awk -v OFS='\t' -F'\t' '{print $1"|"$2, $5}' $outputDir/$outPrefix.pileup.txt |sort -t$'\t') \
|awk -v OFS='\t' -F'\t' '{split($1,a,"|"); print(a[1],a[2], $2,$3,$4,$5,$6,$7, $8,$9,$10)}' \
|awk -v OFS='\t' -F'\t' '{if($11=="") $11 = "NA"; print $0}' > $outputDir/$outPrefix.depth.txt

join -t$'\t' -a 1 -1 1 -2 1 \
<(awk -v OFS='\t' -F'\t' '{print $1"|"$2"|"$3"|"$4, $5,$6,$7,$8,$9,$10,$11}' $outputDir/$outPrefix.depth.txt |sort -t$'\t') \
<(awk -v OFS='\t' -F'\t' '{print $1"|"$2"|"$3"|"$4, $6,$7,$8,$9,$10}' $outputDir/$outPrefix.pileup.txt |sort -t$'\t') \
|awk -v OFS='\t' -F'\t' '{split($1,a,"|"); print(a[1],a[2],a[3],a[4], $2,$3,$4,$5,$6,$7, $8,$9,$10,$11,$12,$13)}' |sort |uniq > $outputDir/$outPrefix.final.txt

sed -i '1i#chr\tpos\tref\talt\tID\tgene\tcDNA\tcodon\ttype\tstrand\tdepth\tref+\tref-\talt+\talt-\tstrandBiasPval'  $outputDir/$outPrefix.final.txt

/opt/fastqc/FastQC/fastqc -o $outputDir -b $bam

rm $outputDir/$outPrefix.depth.txt 

module load musket
module load freebayes
module load samtools

fastq=$1
outDir=$2
reference=$3

fileName=${fastq##*/}
sample=${fileName%%.*}
referenceFile=${reference##*/}
referenceName=${referenceFile%%.*}



#####trimming fastq#####
java -jar /opt/trimmomatic/Trimmomatic-0.33/trimmomatic-0.33.jar \
SE -threads 22 \
$fastq $outDir/$sample.trimmed.fastq \
ILLUMINACLIP:/opt/trimmomatic/Trimmomatic-0.33/NexteraPE-PE.fasta:2:30:10:2:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#####Error correction#####
musket -p 22 $outDir/$sample.trimmed.fastq -o $outDir/$sample.trimmed.ec.fastq
	

	#####SRST2#####
#srst2 --input_se $j"_ec.fastq.gz" \
#--output $j \
#--log \
#--mlst_db /data/SRST2/GAS_MLST_DB/Streptococcus_pyogenes.fasta \
#--mlst_definitions /data/SRST2/GAS_MLST_DB/spyogenes.txt \
#--mlst_delimiter "-" --gene_db /data/SRST2/GAS_emm-plus_DB/emm-plus.fasta

echo "###########aligning"
#####SMALT#####
reference_index=${reference%.fasta}
smalt map -x -n 22 -o $outDir/$sample.sam $reference_index $outDir/$sample.trimmed.ec.fastq
samtools view -bS $outDir/$sample.sam > $outDir/$sample.bam
samtools sort $outDir/$sample.bam $outDir/$sample.sorted
samtools index $outDir/$sample.sorted.bam
samtools flagstat $outDir/$sample.sorted.bam > $outDir/$sample.sorted.bam.flagstat

####Necessar??#### extract unaligned reads from .bam ##############
	##bamtools filter -isMapped false -in $k"-"$j".bam" -out $k"-"$j"_unalnd.bam";
	##bamtools convert -format fastq -in $k"-"$j"_unalnd.bam" -out $k"-"$j"_unalnd.fastq";


#####FreeBayes#####

freebayes -p 1 -= -F 0.7 -f $reference $outDir/$sample.sorted.bam -v $outDir/$sample.vcf
	

#####vcf filtering#####	

/opt/vcflib/vcflib/bin/vcffilter -f "QUAL > 30" $outDir/$sample.vcf |\
/opt/vcflib/vcflib/bin/vcffilter -f "DP > 9" |\
/opt/vcflib/vcflib/bin/vcfallelicprimitives > $outDir/$sample.filtered.vcf

##Need alternatives for vcfstats, as the code is broken	
##/opt/vcflib/vcflib/bin/vcfstats $outDir/$sample.vcf > $outDir/$sample.vcf.stats
##/opt/vcflib/vcflib/bin/vcfstats $outDir/$sample.filtered.vcf > $outDir/$sample.filtered.vcf.stats
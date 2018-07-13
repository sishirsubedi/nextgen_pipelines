
if test $# -gt 0
	then
	while getopts :f:n:b:o: opt
	do
	case $opt in
	f)
		fastq=$OPTARG
		;;
	n)
		readNumber=$OPTARG
		;;
	b)
		database=$OPTARG
		;;
	o)
		outDir=$OPTARG
		;;
	:)
		echo "Option -$OPTARG requires an argument."
		;;
	\?)
		echo "Invalid option: -$OPTARG"
	esac
	done
	shift $((OPTIND-1))


echo "###################################################"
echo "processing " $fastq
#get current time
currentTime=$(date)
currentTimeS1=$(date +%s)
echo "Started at " $currentTime
fileName=${fastq##*/}
sample=${fileName%%.fastq*}
fastqLine=$((readNumber*4))



###FASTQ file extracting and trimming
if [ "$readNumber" != "0" ]
then
gunzip -c $fastq |head -n $fastqLine > $outDir/$sample.trimmed$readNumber.fastq
elif [ "$readNumber" == "0" ]
then
gunzip -c $fastq > $outDir/$sample.trimmed$readNumber.fastq
fi


### FASTQ to FASTA Converstion

echo "Fastq to fasta conversion"
perl /home/MDXLAB/bin/fq2fa.pl fq2fa $outDir/$sample.trimmed$readNumber.fastq > $outDir/$sample.trimmed$readNumber.fasta

### Blast
echo "Blasting"
if [ "$database" != "HN" ]
then
echo "blast against the $database database"
/home/MDXLAB/bin/blastn -query $outDir/$sample.trimmed$readNumber.fasta -db $database -negative_gilist /home/MDXLAB/bin/sequence.gi.txt -num_threads 5 -num_descriptions 1 -num_alignments 1 > $outDir/$sample.blast$readNumber.txt
else
echo "blast against the flu database"
/home/MDXLAB/bin/blastn -query $outDir/$sample.trimmed$readNumber.fasta -db $database -num_threads 5 -num_descriptions 1 -num_alignments 1 > $outDir/$sample.blast$readNumber.txt
fi
### check if Blast process completed
if grep -q "Gap Penalties" $outDir/$sample.blast$readNumber.txt
then
echo "Blast completed"
else
echo "Warning: Blast did not complete correctly"
fi

### Contig counting
echo "Contig Counting"
if [ "$database" != "HN" ]
then
python /home/MDXLAB/bin/contigcount.py -kf 2,3 -exclude /home/MDXLAB/bin/exclude.txt $outDir/$sample.blast$readNumber.txt > $outDir/$sample.count$readNumber.txt
python /home/MDXLAB/bin/contigcount.py -kf 2,3,4 -exclude /home/MDXLAB/bin/exclude.txt $outDir/$sample.blast$readNumber.txt > $outDir/$sample.count$readNumber.extend.txt
else
python /home/MDXLAB/bin/contigcount.py -kf 1,2 -exclude /home/MDXLAB/bin/exclude.txt $outDir/$sample.blast$readNumber.txt > $outDir/$sample.count$readNumber.txt
python /home/MDXLAB/bin/contigcount.py -kf 1,2,3 -exclude /home/MDXLAB/bin/exclude.txt $outDir/$sample.blast$readNumber.txt > $outDir/$sample.count$readNumber.extend.txt
fi
### Complete
currentTime=$(date)
currentTimeS2=$(date +%s)
echo "Completed at " $currentTime

#calculate elapsed time
elapsed=$((currentTimeS2-currentTimeS1))
hours=$((elapsed/3600))
hour_remain=$((elapsed%3600))
minutes=$((hour_remain/60))
seconds=$((hour_remain%60))

echo "Elapsed Time: " $hours " hours" $minutes " minutes" $seconds " seconds" 


else
	echo "Take an input zipped fastq file, and perform blast and contig count."
	echo "Output files will be in the working directory."
	echo "Usage: sh blastn.sh -b [database name] 
	-f [input zipped fastq file] 
	-n [number of reads to blast, input 0 for full blast]
	-o [output directory]"
fi


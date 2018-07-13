
if test $# -gt 0
	then
	while getopts :d:n:b:z: opt
	do
	case $opt in
	d)
		dir=$OPTARG
		;;
	n)
		readNumber=$OPTARG
		;;
	b)
		database=$OPTARG
		;;
	z)
		zipped=$OPTARG
		;;
	:)
		echo "Option -$OPTARG requires an argument."
		;;
	\?)
		echo "Invalid option: -$OPTARG"
	esac
	done
	shift $((OPTIND-1))

if [ $zipped == "yes" ]
then
fileList=$dir/*.fastq.gz
else
fileList=$dir/*.fastq
fi

for file in $fileList
do
echo "###################################################"
echo "processing " $file
#get current time
currentTime=$(date)
currentTimeS1=$(date +%s)
echo "Started at " $currentTime
fileName=${file##*/}
sample=${fileName%%.fastq*}
fastqLine=$((readNumber*4))

###FASTQ file trimming
if [ $readNumber != "0" ] && [ $zipped == "yes" ]
then
gunzip -c $file |head -n $fastqLine > $sample.trimmed$readNumber.fastq
elif [ $readNumber == "0" ] && [ $zipped == "yes" ]
then
gunzip -c $file > $sample.trimmed$readNumber.fastq
elif [ $readNumber != "0" ] && [ $zipped == "no" ]
then
head -n $fastqLine $file > $sample.trimmed$readNumber.fastq
elif [ $readNumber == "0" ] && [ $zipped == "no" ]
then
ln -s $file > $sample.trimmed$readNumber.fastq
fi


### FASTQ to FASTA Converstion

echo "Fastq to fasta conversion"
perl /home/MDXLAB/bin/fq2fa.pl fq2fa $sample.trimmed$readNumber.fastq > $sample.trimmed$readNumber.fasta

### Blast
echo "Blasting"
blastn -query $sample.trimmed$readNumber.fasta -db $database -negative_gilist /home/MDXLAB/bin/sequence.gi.txt -num_threads 6 -num_descriptions 1 -num_alignments 1 > $sample.blast$readNumber.txt

### check if Blast process completed
if grep -q "Gap Penalties" $sample.blast$readNumber.txt
then
echo "Blast completed"
else
echo "Warning: Blast did not complete correctly"
fi

### Contig counting
echo "Contig Counting"
python /home/MDXLAB/bin/contigcount.py -kf 2,3 -exclude /home/MDXLAB/bin/exclude.txt $sample.blast$readNumber.txt > $sample.count$readNumber.txt

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
done

else
	echo "Take an input directory containing zipped fastq files, and perform blast and contig count for all files in the directory."
	echo "Output files will be in the working directory."
	echo "Usage: sh blastn.sh -b [database name] 
	-d [input directory containing zipped fastq files] 
	-n [number of reads to blast, input 0 for full blast] 
	-z [whether the input files are zipped, yes or no]"
fi
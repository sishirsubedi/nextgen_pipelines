if test $# -gt 0
	then
	while getopts :p:q:o: opt
	do
	case $opt in
	p)
		fastq1=$OPTARG
		;;
	q)
		fastq2=$OPTARG
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


fileName=${fastq1##*/}
sample=${fileName%%_*}

###align###
#sh /home/pipelines/master/shell/exomeAlign.sh $fastq1 $fastq2 $outDir
###gatkScan###
sh /home/pipelines/master/shell/gatkScan.sh "$outDir"/"$sample".sort.bam $outDir


fi

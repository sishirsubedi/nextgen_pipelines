###Takes a 50 gene excel file and load into table molseq.gene50#####

##make working directory in /home/scratch


if test $# -gt 0
	then
	while getopts :r:b: opt
	do
	case $opt in
	r)
		report_no=$OPTARG
		;;
	b)
		barcode=$OPTARG
		;;
	:)
		echo "Option -$OPTARG requires an argument."
		;;
	\?)
		echo "Invalid option: -$OPTARG"
	esac
	done
	shift $((OPTIND-1))	

#make dir in /home/scratch

if [ ! -d /home/scratch/$report_no-$barcode ]
then
mkdir /home/scratch/$report_no-$barcode
fi
	
#excel ---> dbsnp file
#locate excel file

infile=$(ls /home/ionadmin/archivedReports/*$report_no/plugin_out/variantCaller_out/IonXpress_$barcode/alleles.xls)
echo "converting to db load file"
python /var/www/html/hotspot/gene50dbLoadweb.py gene50excelLoad \
-I $infile \
-o /home/scratch/$report_no-$barcode/alleles.db.txt \
-r $report_no \
-b $barcode

#prep amplicon coverage loading
echo "Preparing the amplicon coverage file to be loaded"

coverageFile=$(ls /home/ionadmin/archivedReports/*$report_no/plugin_out/coverageAnalysis_out/IonXpress_$barcode/*.amplicon.cov.xls)

if [ -z $coverageFile ]
then
echo "coverage file not found"
else
echo " loading coverage file into database"
awk -v OFS='\t' '{print $1, $2,$3, $4, $5, $10, "'$report_no'", "'$barcode'"}'  $coverageFile |tail -n +2 > /home/scratch/$report_no-$barcode/coverage.db.txt
mysql -e "load data local infile '"/home/scratch/$report_no-$barcode/coverage.db.txt"' into table molseq.amplicon_coverage \
(chr, start, end, amplicon, gene, totalReads, report_no, barcode)" -u webuser --password=molSeq3127
fi
#load into database
echo "database loading"
mysql -e "load data local infile '"/home/scratch/$report_no-$barcode/alleles.db.txt"' into table molseq.gene50 \
(chr, pos, cosmic, ref, alt, genotype, depth, altCount, frequency, report_no, barcode)" -u webuser --password=molSeq3127

echo "Finished"


else
	echo "Usage: sh gene50load.sh -r [run number] -b [barcode] "
fi



###Takes a 50 gene excel file and load into table molseq.gene50#####

if test $# -gt 0
	then
	while getopts :r:v:c:b: opt
	do
	case $opt in
	r)
		report_no=$OPTARG
		;;
	v)
		variantCaller=$OPTARG
		;;
	c)
		coverageID=$OPTARG
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
	

#excel ---> dbsnp file
#locate excel file
infile=$(ls /home/proton/*$report_no/plugin_out/variantCaller_out.$variantCaller/IonXpress_$barcode/alleles.xls)
#mkdir for the dbload files
report_sur=${infile#/home/proton/}
report=${report_sur%/plugin_out/variantCaller_out.$variantCaller/IonXpress_$barcode/alleles.xls}
echo "making directories in /home/proton/dbload"
if [ ! -d "/home/proton/dbload/$report" ]
then
mkdir /home/proton/dbload/$report
fi

if [ ! -d "/home/proton/dbload/$report/IonXpress_$barcode" ]
then
mkdir /home/proton/dbload/$report/IonXpress_$barcode
fi

echo "converting to db load file"
python /home/niyunyun/code/python/gene50dbLoadweb.py gene50excelLoad \
-I $infile \
-o /home/proton/dbload/$report/IonXpress_$barcode/alleles.db.txt \
-r $report_no \
-b $barcode


#load into database
echo "database loading"
mysql -e "load data local infile '"/home/proton/dbload/$report/IonXpress_$barcode/alleles.db.txt"' into table molseq.gene50_proton \
(chr, pos, cosmic, ref, alt, genotype, depth, altCount, frequency, report_no, barcode)" -u niyunyun --password=molSeq3127

#prep amplicon coverage loading
echo "Preparing the amplicon coverage file to be loaded"

coverageFile=$(ls /home/proton/*$report_no/plugin_out/coverageAnalysis_out.$coverageID/IonXpress_$barcode/*.amplicon.cov.xls)

if [ -z $coverageFile ]
then
echo "coverage file not found"
else
echo " loading coverage file into database"
awk -v OFS='\t' '{print $1, $2,$3, $4, $5, $10, "'$report_no'", "'$barcode'"}'  $coverageFile |tail -n +2 > /home/proton/dbload/$report/IonXpress_$barcode/coverage.db.txt
mysql -e "load data local infile '"/home/proton/dbload/$report/IonXpress_$barcode/coverage.db.txt"' into table molseq.amplicon_coverage_proton \
(chr, start, end, amplicon, gene, totalReads, report_no, barcode)" -u webuser --password=molSeq3127
fi


echo "Finished"


else
	echo "Usage: sh gene50load_proton.sh -r [report number] -b [barcode] "
fi



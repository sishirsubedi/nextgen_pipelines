###Takes a 50 gene excel file and load into table molseq.gene50#####

if test $# -gt 0
	then
	while getopts :r:b:c:v: opt
	do
	case $opt in
	r)
		report_no=$OPTARG
		;;
	b)
		barcode=$OPTARG
		;;
	v)
		callerID=$OPTARG
		;;
	c)
		coverageID=$OPTARG
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

#mkdir for the dbload files
if [ ! -z $callerID ]
then
infile=$(ls /home/ionadmin/archivedReports/*$report_no/plugin_out/variantCaller_out.$callerID/IonXpress_$barcode/alleles.xls)
report_sur=${infile#/home/ionadmin/archivedReports/}
report=${report_sur%/plugin_out/variantCaller_out.$callerID/IonXpress_$barcode/alleles.xls}
else
infile=$(ls /home/ionadmin/archivedReports/*$report_no/plugin_out/variantCaller_out/IonXpress_$barcode/alleles.xls)
report_sur=${infile#/home/ionadmin/archivedReports/}
report=${report_sur%/plugin_out/variantCaller_out/IonXpress_$barcode/alleles.xls}
fi

echo "making directories in /home/ionadmin/dbload"
if [ ! -d "/home/ionadmin/dbload/$report" ]
then
mkdir /home/ionadmin/dbload/$report
fi

if [ ! -d "/home/ionadmin/dbload/$report/IonXpress_$barcode" ]
then
mkdir /home/ionadmin/dbload/$report/IonXpress_$barcode
fi

echo "converting to db load file"
if [ ! -z $callerID ]
then
python /home/niyunyun/code/python/gene50dbLoadweb.py gene50excelLoad \
-I $infile \
-o /home/ionadmin/dbload/$report/IonXpress_$barcode/alleles.db.$callerID.txt \
-r $report_no \
-b $barcode \
-c $callerID
echo "database loading"
mysql -e "load data local infile '"/home/ionadmin/dbload/$report/IonXpress_$barcode/alleles.db.$callerID.txt"' into table molseq.gene50 \
(chr, pos, cosmic, ref, alt, genotype, depth, altCount, frequency, report_no, barcode,callerID)" -u niyunyun --password=molSeq3127 
else
python /home/niyunyun/code/python/gene50dbLoadweb.py gene50excelLoad \
-I $infile \
-o /home/ionadmin/dbload/$report/IonXpress_$barcode/alleles.db.txt \
-r $report_no \
-b $barcode \
-c "None"
echo "database loading"
mysql -e "load data local infile '"/home/ionadmin/dbload/$report/IonXpress_$barcode/alleles.db.txt"' into table molseq.gene50 \
(chr, pos, cosmic, ref, alt, genotype, depth, altCount, frequency, report_no, barcode,callerID)" -u niyunyun --password=molSeq3127 
fi



#prep amplicon coverage loading
echo "Preparing the amplicon coverage file to be loaded"



#if [ -z $coverageFile ]
#then
#echo "coverage file not found"
#else
#echo " loading coverage file into database"
if [ ! -z $coverageID ]
then
coverageFile=$(ls /home/ionadmin/archivedReports/*$report_no/plugin_out/coverageAnalysis_out.$coverageID/IonXpress_$barcode/*.amplicon.cov.xls)
awk -v OFS='\t' '{print $1, $2,$3, $4, $5, $10, "'$report_no'", "'$barcode'", "'$coverageID'"}'  $coverageFile |tail -n +2 > /home/ionadmin/dbload/$report/IonXpress_$barcode/coverage.db$coverageID.txt
else
coverageFile=$(ls /home/ionadmin/archivedReports/*$report_no/plugin_out/coverageAnalysis_out/IonXpress_$barcode/*.amplicon.cov.xls)
awk -v OFS='\t' '{print $1, $2,$3, $4, $5, $10, "'$report_no'", "'$barcode'", "None"}'  $coverageFile |tail -n +2 > /home/ionadmin/dbload/$report/IonXpress_$barcode/coverage.db$coverageID.txt
fi

mysql -e "load data local infile '"/home/ionadmin/dbload/$report/IonXpress_$barcode/coverage.db$coverageID.txt"' into table molseq.amplicon_coverage \
(chr, start, end, amplicon, gene, totalReads, report_no, barcode, coverageID)" -u webuser --password=molSeq3127



echo "Finished"


else
	echo "Usage: sh gene50load.sh -r [report number] -b [barcode] -v [variant caller ID] -c [coverage ID] "
fi



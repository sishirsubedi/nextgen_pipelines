###generate 50 gene report#####

if test $# -gt 0
	then
	while getopts :r:b:v:c: opt
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
	
#load results into database
if [ ! -z $callerID ]
then
sh /home/niyunyun/code/shell/gene50load.sh -r $report_no -b $barcode -v $callerID -c $coverageID
else
sh /home/niyunyun/code/shell/gene50load.sh -r $report_no -b $barcode
fi

#Query and generate report#####
if [ ! -z $callerID ]
then
bash ~/code/shell/gene50dbQuery.sh -r $report_no -b $barcode -v $callerID -c $coverageID
else
bash ~/code/shell/gene50dbQuery.sh -r $report_no -b $barcode -v "None" -c "None"
fi

##generate the human readable report###
python /home/niyunyun/code/python/gene50ParseFilter.py \
-I /home/gene50Report/$report_no-$barcode.out.filtered.txt \
-o /home/gene50Report/$report_no-$barcode.out.filtered.report.txt \
-r $report_no -b $barcode \
-c /home/gene50Report/$report_no-$barcode.out.amplicon.txt
else
	echo "Usage: sh gene50report.sh -r [report number] -b [barcode] -v [variant caller ID] -c [coverage analysis ID] "
fi

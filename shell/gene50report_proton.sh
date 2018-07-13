###generate 50 gene report#####

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
	
#load results into database
sh /home/niyunyun/code/shell/gene50load_proton.sh -r $report_no -b $barcode -v $variantCaller -c $coverageID

#Query and generate report#####
bash ~/code/shell/gene50dbQuery_proton.sh -r $report_no -b $barcode

else
	echo "Usage: sh gene50report.sh -r [report number] -v [variant caller ID] -c [coverage analysis ID] -b [barcode] "
fi

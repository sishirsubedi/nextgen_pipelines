###Takes a vcf file and load into table molseq.breast_her#####

if test $# -gt 0
	then
	while getopts :v:s:o: opt
	do
	case $opt in
	v)
		vcf=$OPTARG
		;;
	s)
		sample=$OPTARG
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
	

#vcf ---> dbsnp file

echo "converting to db load file"
prefix=${vcf%.vcf}
python ~/code/python/vcfdbLoad.py brca_her_load \
-I $vcf -o $outDir/$prefix.filter.txt -r /doc/ref/clinvar_20150305.vcf -s $sample

#run snpEff
echo "running snpEff"
java -Xmx4g -jar /opt/software/snpEff-4.1b/snpEff.jar GRCh37.75 $outDir/$prefix.filter.txt > $outDir/$prefix.filter.eff.txt

#parse snpEff to dbload file
echo "producing dbload file"
python ~/code/python/vcfdbLoad.py eff_to_db -I $outDir/$prefix.filter.eff.txt -o $outDir/$prefix.dbload.txt

#load into database
echo "database loading"
mysql -e "load data local infile '"$outDir/$prefix.dbload.txt"' into table molseq.breast_her ignore 2 lines" -u niyunyun --password=molSeq3127




else
	echo "Usage: sh breast_her_load.sh -v [input vcf file] -s [sample name] -o [out dir]"
fi



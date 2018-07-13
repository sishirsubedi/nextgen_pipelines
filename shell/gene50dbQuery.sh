###query 50gene database#####

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
	



#Query database
echo "database query"
mysql -e \
"select distinct t1.*, \
t2.cDNA, t2.protein, t2.coding, t2.gene, \
t3.id as dbsnp, t3.property as property_dbSNP, 
REPLACE(REPLACE(t4.property, '\r', ''), '\n', ''), REPLACE(REPLACE(t4.note, '\r', ''), '\n', ''), \
t6.occurance_count \
into outfile '"/home/scratch/$report_no-$barcode.out.txt"' \
fields terminated by '\t' \
lines terminated by '\n' \
from gene50 as t1 \
left join cosmic as t2 \
on t1.cosmic = t2.id \
left join dbsnp as t3 \
on t1.chr = t3.chr and \
t1.pos = t3.pos and \
t1.ref = t3.ref and \
t1.alt = t3.alt \
left join \
gene50_curation as t4 \
on t1.cosmic = t4.cosmic \
left join \
(select cosmic, count(*) as occurance_count from gene50 \
where (genotype = 'Homozygous' or genotype = 'Heterozygous') and cosmic != '---' \
group by cosmic) as t6 \
on t1.cosmic = t6.cosmic \
where report_no = '"$report_no"' and barcode = '"$barcode"' and callerID = '"$callerID"' \
order by t1.cosmic desc" \
-u niyunyun -p molseq --password=molSeq3127 

#mysql -e \
#"select distinct t1.*, \
#t2.cDNA, t2.protein, t2.coding, t2.gene, \
#t3.id as dbsnp, t3.property as property_dbSNP,
#REPLACE(REPLACE(t4.property, '\r', ''), '\n', ''), REPLACE(REPLACE(t4.note, '\r', ''), '\n', ''),  \
#t6.occurance_count \
#into outfile '"/home/scratch/$report_no-$barcode.out.filtered.txt"' \
#fields terminated by '\t' \
#lines terminated by '\n' \
#from gene50 as t1 \
#left join cosmic as t2 \
#on t1.cosmic = t2.id \
#left join dbsnp as t3 \
#on t1.chr = t3.chr and \
#t1.pos = t3.pos and \
#t1.ref = t3.ref and \
#t1.alt = t3.alt \
#left join \
#gene50_curation as t4 \
#on t1.cosmic = t4.cosmic \
#left join \
#(select cosmic, count(*) as occurance_count from gene50 \
#where (genotype = 'Homozygous' or genotype = 'Heterozygous') and cosmic != '---' \
#group by cosmic) as t6 \
#on t1.cosmic = t6.cosmic \
#where report_no = '"$report_no"' and barcode = '"$barcode"' and \
#t1.cosmic != '---' and frequency > 9.0 and coding != 'synonymous' and coding != 'non-coding' and t1.genotype != 'No Call' \
#order by t1.cosmic desc" \
#-u niyunyun -p molseq --password=molSeq3127 


#move files into the report folder
#mv /home/scratch/$report_no-$barcode.out.filtered.txt /home/gene50Report
mv /home/scratch/$report_no-$barcode.out.txt /home/gene50Report

#filter the output
awk -F'\t' '{
if(($4 != "---" &&
$10 > 9.0 &&
$17 != "synonymous" &&
$17 != "non-coding" &&
$7 != "No Call") ||
($10 > 9.0 &&
$17 != "synonymous" &&
$17 != "non-coding" &&
$7 != "No Call" &&
length($5) != length($6))) print}' /home/gene50Report/$report_no-$barcode.out.txt > /home/gene50Report/$report_no-$barcode.out.filtered.txt


#combine raw and filtered files
cat <(sed -e "\$aFiltered\n" /home/gene50Report/$report_no-$barcode.out.txt) /home/gene50Report/$report_no-$barcode.out.filtered.txt > /home/gene50Report/$report_no-$barcode.out.final.txt

#query amplicon database
mysql -e \
"select t1.*, t2.hotspotID \
into outfile '"/home/scratch/$report_no-$barcode.out.amplicon.txt"' \
fields terminated by '\t' \
lines terminated by '\n' \
from amplicon_coverage as t1 \
left join gene50_amplicon as t2 \
on t1.amplicon = t2.name \
where report_no = '"$report_no"' and barcode = '"$barcode"' and coverageID = '"$coverageID"' and \
totalReads < 100" \
-u niyunyun -p molseq --password=molSeq3127 

#copy amplicon file to report dir
mv /home/scratch/$report_no-$barcode.out.amplicon.txt /home/gene50Report
#add header to amplicon file
sed -i '1ichr\tstart\tend\tamplicon\tgene\ttotalReads\treport_no\tbarcode\tcovAnalysisID\tmutations' /home/gene50Report/$report_no-$barcode.out.amplicon.txt

#add header line to output files

sed -i '1iuid\tchr\tpos\tcosmic\tref\talt\tgenotype\tdepth\taltCount\tfreqeuncy\treport_no\tbarcode\tcallerID\treportStatus\tcDNA\tprotein\tchangeType\tgene\tdbSNP\tdbSNPAnnotation\tproperty\tnote\toccurance' \
/home/gene50Report/$report_no-$barcode.out.txt



#combine raw and filtered files
cat \
<(sed -e "\$aFiltered\n" /home/gene50Report/$report_no-$barcode.out.txt) \
/home/gene50Report/$report_no-$barcode.out.filtered.txt \
>/home/gene50Report/$report_no-$barcode.out.combineFiltered.txt
cat \
<(sed -e "\$aInsufficiently Covered Amplicons\n" /home/gene50Report/$report_no-$barcode.out.combineFiltered.txt) \
/home/gene50Report/$report_no-$barcode.out.amplicon.txt \
> /home/gene50Report/$report_no-$barcode.out.combine.txt

#remove line mysql line break in output file
awk '{ sub(/\.\\/,""); print }' /home/gene50Report/$report_no-$barcode.out.combine.txt > /home/gene50Report/$report_no-$barcode.out.final.txt

#remove intermediate files

rm /home/gene50Report/$report_no-$barcode.out.combineFiltered.txt
rm /home/gene50Report/$report_no-$barcode.out.combine.txt

echo "Finished"


else
	echo "Usage: sh gene50dbQuery.sh -r [report number] -b [barcode] "
fi

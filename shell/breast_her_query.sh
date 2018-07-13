###database query and output####



if test $# -gt 0
	then
	while getopts :p:s:o: opt
	do
	case $opt in
	p)
		pathogenic=$OPTARG
		;;
	s)
		sample=$OPTARG
		;;
	o)
		output=$OPTARG
		;;
	:)
		echo "Option -$OPTARG requires an argument."
		;;
	\?)
		echo "Invalid option: -$OPTARG"
	esac
	done
	shift $((OPTIND-1))

if [ $pathogenic == "y" ]

then
mysql_cmd="
select t1.*,t2.gene, t2.clnsig, t2.clnsig_db, t2.clnsig_disease, t2.clnallele, t2.caf, t2.common, t2.curated_clnsig  
from molseq.breast_her t1
left join molseq.clinvar t2
on (t1.chr = t2.chr and t1.pos = t2.pos and
t1.dbSNP = t2.dbSNP and t1.ref = t2.ref and
t1.alt = t2.alt)
where sample = '"$sample"' and frequency != '0.00' and (clnsig like '%pathogenic%' or clnsig like '%likely_pathogenic%' or clnsig like '%uncertain_significance%')
into outfile '"output/$output"'"

else
mysql_cmd="
select t1.*,t2.gene, t2.clnsig, t2.clnsig_db, t2.clnsig_disease, t2.clnallele, t2.caf, t2.common, t2.curated_clnsig  
from molseq.breast_her t1
left join molseq.clinvar t2
on (t1.chr = t2.chr and t1.pos = t2.pos and
t1.dbSNP = t2.dbSNP and t1.ref = t2.ref and
t1.alt = t2.alt)
where sample = '"$sample"' and frequency != '0.00' 
into outfile '"output/$output"'"

fi

#query database
mysql -e "$mysql_cmd" -u niyunyun --password=molSeq3127

#move output to current directory
mv /var/lib/mysql/output/$output ./

#add header to output
sed -i '1i@chr\tpos\tdbSNP\tref\talt\tfrequency\ttype\tgenotype\tmulti\tfixed\tbase_corrected\tsample\tgeneID\tchangeType\tcDNAchange\tproteinChange\tgene\tclnsig\tclnsig_db\tclnsig_disease\tclnallele\tcaf\tcommon\tcurated_clnsig' $output

else
	echo "Usage: sh breast_her_query.sh -s [sample] -o [output file, name only, will be saved in current dir] -p [y: for pathogenic and undetermined results only; n: for all results]"
fi
#!/bin/bash
if test $# -gt 0
	then
	while getopts :s:i:o: opt
	do
	case $opt in
	s)
		snp=$OPTARG
		;;
	i)
		indel=$OPTARG
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
	
	
	fileName=${snp##*/}
	sample=${fileName%%.*}
	echo "Post processing Varscan results for $sample"	
	#####combine snp and indel calls#####
	echo "combine snp and indel calls"

cat <(awk -v OFS='\t' '{print $1,$2,$3,$19, $5, $6, $10, $11}' $snp |tail -n +2) \
<(awk -v OFS='\t' '{
alt=substr($19, 2)
if($19~/^-/) 
print $1,$2,$3alt, $3, $5, $6, $10, $11
else if($19~/^+/)
print $1, $2, $3, $3alt, $5, $6, $10, $11}' $indel) |\
sort -k1,1 -k2,2n > $outDir/$sample.comb.txt
	
	####convert combined txt file into vcf file####
	echo "converting to vcf"
	python /home/niyunyun/code/python/varScan2vcf.py -I $outDir/$sample.comb.txt -o $outDir/$sample.comb.vcf
	
else
	echo "Usage: sh postVarScan.sh -s [VarScan SNP result] -i [VarScan indel result] -o [output directory]"
fi

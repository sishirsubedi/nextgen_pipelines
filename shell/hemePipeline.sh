
file=$1
file_name=${file##*/}
dir=${file%/*}
run=${dir##*/}
#echo $dir
if [ $file_name == "CompletedJobInfo.xml" ]
then
#run filter vcf
sleep 15m
echo "Start post processing ..." >> $dir/Data/Intensities/BaseCalls/Alignment/postProcess.log
sampleNum=$(ls $dir/Data/Intensities/BaseCalls/Alignment/*.vcf |wc -l)
echo "Total samples to process: "$sampleNum >> $dir/Data/Intensities/BaseCalls/Alignment/postProcess.log
if [[ $sampleNum -ge 1 ]]
then
for file in $dir/Data/Intensities/BaseCalls/Alignment/*.vcf
do
fileName=${file##*/}
sample=${fileName%.vcf}
sampleShort=${sample%%_*}
sampleLower=${sampleShort,,}
echo "Start processing $sample" >> $dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log 2>&1
grep "#CHROM" $file > $dir/Data/Intensities/BaseCalls/Alignment/$sample.header
echo "filtering out excluded exons" >> $dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log 2>&1
/opt/software/bedtools-2.17.0/bin/bedtools intersect -u -a $file -b /home/doc/ref/Heme/trusight-myeloid-amplicon-track.excluded.bed > $dir/Data/Intensities/BaseCalls/Alignment/$sample.amp.vcf 2>>$dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log
cat $dir/Data/Intensities/BaseCalls/Alignment/$sample.header $dir/Data/Intensities/BaseCalls/Alignment/$sample.amp.vcf > $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.vcf 2>>$dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log
##read depth >= 100
echo "filtering out read depth < 100" >> $dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log 2>&1
python /home/niyunyun/code/python/filterVcf.py \
-I $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.vcf \
-o $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filterDP.vcf \
-i -f DP -p ">=" -v "100" >> $dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log 2>&1
##alt freq >= 10
if [ $sampleLower = "horizondna" ]
then
echo "filtering out alt freq < 1% for horizon" >> $dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log 2>&1
python /home/niyunyun/code/python/filterVcf.py \
-I $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filterDP.vcf \
-o $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filterDP.filterFreq.vcf \
-n 10 -f VF -p ">=" -v "0.01" >> $dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log 2>&1
else
echo "filtering out alt freq < 10% for others" >> $dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log 2>&1
python /home/niyunyun/code/python/filterVcf.py \
-I $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filterDP.vcf \
-o $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filterDP.filterFreq.vcf \
-n 10 -f VF -p ">=" -v "0.1" >> $dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log 2>&1
fi
##filter out not PASS
echo "filtering out NOT pass" >> $dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log 2>&1
awk '$7=="PASS"' $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filterDP.filterFreq.vcf > $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filterDP.filterFreq.pass.vcf 2>>$dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log
cat $dir/Data/Intensities/BaseCalls/Alignment/$sample.header $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filterDP.filterFreq.pass.vcf > $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filter.vcf 2>>$dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log
#run VEP
echo "Running VEP" >> $dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log 2>&1
/home/niyunyun/perl/bin/perl \
/opt/vep/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \
-i $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filter.vcf \
-o $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filter.vep.vcf \
--offline \
--dir_cache /opt/vep/ensembl-tools-release-83/cache/ \
--sift p \
--polyphen p \
--hgvs \
--symbol \
--vcf \
--pubmed \
--fasta /home/doc/ref/ref_genome/Homo_sapiens.GRCh37.73.dna_sm.primary_assembly.fa \
--pick_allele \
--individual all \
--check_alleles \
--force_overwrite \
--gmaf \
--maf_1kg >> $dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log 2>&1

#parse vep result
echo "Parsing VEP" >> $dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log 2>&1
python /home/niyunyun/code/python/parseVEP.py parseHeme \
-I $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filter.vep.vcf \
-o $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filter.vep.parse.txt \
>> $dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log 2>&1

#filter on VEP impacts
echo "filtering on VEP impacts" >> $dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log 2>&1
awk '{if($6=="HIGH" || $6 =="MODERATE") print}' $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filter.vep.parse.txt > $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filter.vep.parse.filterImpact.txt 2>>$dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log

###join with clinvar
echo "Joining with Clinvar" >> $dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log 2>&1
join -t$'\t' -a 1 -1 1 -2 1 \
<(awk -F'\t' -v OFS='\t' '{print($2"|"$3"|"$4"|"$5,$1,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28)}' $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filter.vep.parse.filterImpact.txt |sort -t$'\t' -k1,1) \
<(awk -F'\t' -v OFS='\t' '{print($1"|"$2"|"$4"|"$5, $9,$4,$5,$6, $8)}' /home/doc/ref/Heme/clinvar_20150305.parse.single.txt  |sort -t$'\t' -k1,1) \
|awk -F'\t' -v OFS='\t' '{if($26 == "") {$26 = ""; $27=""; $28 = ""; $29 = ""; $30 = ""} print $0}' >  $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filter.vep.parse.filterImpact.clinVar.txt 2>>$dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log

###join with cosmic
echo "Joining with Cosmic" >> $dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log 2>&1
join -t$'\t' -a 1 -1 1 -2 1 \
<(sort -t$'\t' -k1,1 $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filter.vep.parse.filterImpact.clinVar.txt) \
<(awk -F'\t' -v OFS='\t' '{print($1"|"$2"|"$4"|"$5, $3)}' /home/doc/ref/Heme/CosmicCodingMutsV72.inHeme.withChr.vcf | sort -t$'\t' -k1,1) \
|awk -F'\t' -v OFS='\t' '{split($1, a, "|"); print($2, a[1], a[2],a[3],a[4],$3,$4,$5, $6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28, $29,$30, $31, $32, $33, $34)}' \
> $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filter.vep.parse.filterImpact.clinvar.cosmic.txt 2>>$dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log

sed -i '1iGene_0\tchr\tpos\tref\talt\tClassification_23\tType_5\tGenotype_6\tQuality_9\tAlt Variant Freq_13\tRead Depth_14\tAlt Read Depth_15\tConsequence_29\tSift_40\tPolyPhen_41\tHGVSc_43\tHGVSp_44\tdbSNP ID_45\tAncestral Allele_46\tAllele Freq Global Minor_48\tGlobal Minor Allele_49\tAllele Freq Amr_50\tAllele Freq Asn_51\tAllele Freq Af_52\tAllele Freq Eur_53\tCOSMIC ID_58\tCOSMIC Wildtype_59\tCOSMIC Allele_60\tClinVar Accession_70\tClinVar Ref_65\tClinVar Alleles_66\tClinVar Allele Type_67\tClinVar Significance_68\tcosmicV72' $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filter.vep.parse.filterImpact.clinvar.cosmic.txt 2>>$dir/Data/Intensities/BaseCalls/Alignment/$sample.postProcess.log

###make symbolic link in miseq_out folder for the parsed file
ln -s $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filter.vep.parse.filterImpact.clinvar.cosmic.txt /home/vmUser/miseq_out/$run/$sampleShort.newCosmic.txt

rm $dir/Data/Intensities/BaseCalls/Alignment/$sample.amp.vcf
rm $dir/Data/Intensities/BaseCalls/Alignment/$sample.header
rm $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filterDP.vcf
rm $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filterDP.filterFreq.vcf
rm $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filterDP.filterFreq.pass.vcf
rm $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filter.vep.parse.txt
rm $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filter.vep.parse.filterImpact.txt
rm $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.filter.vep.parse.filterImpact.clinVar.txt

done
grep -v -f /home/doc/ref/Heme/excludedAmplicons.txt $dir/Data/Intensities/BaseCalls/Alignment/AmpliconCoverage_M1.tsv > $dir/Data/Intensities/BaseCalls/Alignment/AmpliconCoverage_M1.filtered.tsv
fi





#after the run, remove the dir from the monitor list
sed -i "\:$dir IN_CREATE bash /home/niyunyun/code/shell/hemePipeline.sh \$@/\$#:d" /var/spool/incron/root
fi
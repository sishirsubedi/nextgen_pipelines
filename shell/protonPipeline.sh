if [ $# -eq 0 ]
then
echo "Usage: protonPipeline.sh"
echo "-r runID: three digits, last number of the run"
echo "-c variant caller ID"
echo "-a amplicon coverage analysis ID"
echo "-e environment type ID"
exit
fi

if test $# -gt 0
	then
	while getopts :r:c:a:e: opt
	do
	case $opt in
	r)
		runID=$OPTARG
		;;
	c)
		callerID=$OPTARG
		;;
	a)
		coverageID=$OPTARG
		;;
	e)
		environmentID=$OPTARG
		;;
	:)
		echo "Option -$OPTARG requires an argument."
		;;
	\?)
		echo "Invalid option: -$OPTARG"
	esac
	done
	shift $((OPTIND-1))
fi

if [ ! -d /home/environments/$environmentID/proton/*"$runID"/plugin_out/coverageAnalysis_out."$coverageID"/ ]
then
echo "coverage folder for run ID:$runID; coverage ID:$coverageID not found"
exit
fi

if [ ! -d /home/environments/$environmentID/proton/*"$runID"/plugin_out/variantCaller_out."$callerID"/ ]
then
echo "variant folder for run ID:$runID; variant caller ID:$callerID not found"
exit
fi


variantFolder=$(ls -d /home/environments/$environmentID/proton/*$runID/plugin_out/variantCaller_out."$callerID")
ampliconFolder=$(ls -d /home/environments/$environmentID/proton/*"$runID"/plugin_out/coverageAnalysis_out."$coverageID")
runFolder=$(ls -d /home/environments/$environmentID/proton/*$runID)
runName=${runFolder##*/}

if [ ! -d /home/environments/$environmentID/protonAnalysis/$runName ]
then
mkdir /home/environments/$environmentID/protonAnalysis/$runName
fi
chmod 777 /home/environments/$environmentID/protonAnalysis/$runName

if [ ! -d /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID" ]
then
mkdir /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"
fi
chmod 777 /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"

if [ ! -d /home/environments/$environmentID/protonAnalysis/$runName/coverageAnalysis_out."$coverageID" ]
then
mkdir /home/environments/$environmentID/protonAnalysis/$runName/coverageAnalysis_out."$coverageID"
fi
chmod 777 /home/environments/$environmentID/protonAnalysis/$runName/coverageAnalysis_out."$coverageID"

echo "Runfolder is $runFolder"

##get runDate information##
declare -A months
months=( ["Jan"]="01" ["Feb"]="02" ["Mar"]="03" ["Apr"]="04" ["May"]="05" ["Jun"]="06" ["Jul"]="07" ["Aug"]="08" ["Sep"]="09" ["Oct"]="10" ["Nov"]="11" ["Dec"]="12" )
if [ -f /home/environments/$environmentID/proton/*"$runID"/InitLog.txt ]
then
runDate=$(head -n 1 /home/environments/$environmentID/proton/*"$runID"/InitLog.txt)
year1=$(echo $runDate |cut -d ' ' -f 5)
year=${year1%:}
day=$(echo $runDate |cut -d ' ' -f 3)
monthWord=$(echo $runDate |cut -d ' ' -f 2)
month=${months["$monthWord"]}
date=$year-$month-$day
echo $date > /home/environments/$environmentID/protonAnalysis/$runName/runDate.txt
else
echo "Warning: InitLog.txt file not found, run date will not be entered"
fi

for file in $variantFolder/IonXpress_*/TSVC_variants.vcf
do
echo "Processing $file"
sampleFolder=${file%/*}
sampleName=${sampleFolder##*/}
if [ ! -d /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"/$sampleName ]
then
mkdir /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"/$sampleName
fi
if [ ! -d /home/environments/$environmentID/protonAnalysis/$runName/coverageAnalysis_out."$coverageID"/$sampleName ]
then
mkdir /home/environments/$environmentID/protonAnalysis/$runName/coverageAnalysis_out."$coverageID"/$sampleName
fi
chmod 777 /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"/$sampleName
chmod 777 /home/environments/$environmentID/protonAnalysis/$runName/coverageAnalysis_out."$coverageID"/$sampleName
echo "split multiple alt alleles into different lines"
python /home/pipelines/master/python/splitVcf.py \
-I $file \
-f FAO,FDP,AF \
-o /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"/$sampleName/TSVC_variants.split.vcf

echo "running VEP"
/home/pipelines/master/perl/bin/perl \
/opt/vep/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \
-i /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"/$sampleName/TSVC_variants.split.vcf \
-o /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"/$sampleName/TSVC_variants.split.vep.vcf \
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
--check_alleles \
--force_overwrite \
--gmaf \
--maf_1kg

echo "parse VEP results"
python /home/pipelines/master/python/parseVEP.py \
parseIon \
-I /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"/$sampleName/TSVC_variants.split.vep.vcf \
-o /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"/$sampleName/TSVC_variants.split.vep.parse.txt

echo "parse VEP results for new varView"
python /home/pipelines/master/python/parseVEP.py \
parseIonNewVarView \
-I /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"/$sampleName/TSVC_variants.split.vep.vcf \
-o /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"/$sampleName/TSVC_variants.split.vep.parse.newVarView.txt

echo "filter VEP results"
awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 10 && $11 >=100) print}' /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"/$sampleName/TSVC_variants.split.vep.parse.txt > /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"/$sampleName/TSVC_variants.split.vep.parse.filter.txt
awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 10 && $11 >=100) print}' /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"/$sampleName/TSVC_variants.split.vep.parse.newVarView.txt > /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"/$sampleName/TSVC_variants.split.vep.parse.newVarView.filter.txt

echo "Joining with Clinvar"
join -t$'\t' -a 1 -1 1 -2 1 \
<(awk -F'\t' -v OFS='\t' '{print($3"|"$4"|"$5"|"$6,$1,$2,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24)}' /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"/$sampleName/TSVC_variants.split.vep.parse.filter.txt |sort -t$'\t' -k1,1) \
<(awk -F'\t' -v OFS='\t' '{print($1"|"$2"|"$4"|"$5, $9,$6, $8)}' /home/doc/ref/clinvar/clinvar_20150305.parse.single.collapse.txt  |sort -t$'\t' -k1,1) \
|awk -F'\t' -v OFS='\t' '{if($22 == "") {$22 = ""; $23=""; $24 = "";} print $0}' >  /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"/$sampleName/TSVC_variants.split.vep.parse.filter.clinvar.txt


echo "Joining with Cosmic"
join -t$'\t' -a 1 -1 1 -2 1 \
<(sort -t$'\t' -k1,1 /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"/$sampleName/TSVC_variants.split.vep.parse.filter.clinvar.txt) \
<(awk -F'\t' -v OFS='\t' '{print($1"|"$2"|"$4"|"$5, $3)}' /home/doc/ref/cosmic/CosmicCodingMutsV72.db.sort.collapse.txt | sort -t$'\t' -k1,1) \
|awk -F'\t' -v OFS='\t' '{split($1, a, "|"); print($2, $3, a[1], a[2],a[3],a[4],$4,$5, $6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25)}' \
> /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"/$sampleName/TSVC_variants.final.txt

sed -i '1iGene\texon\tchr\tpos\tref\talt\tClassification\tType\tQuality\tAltVariantFreq\tRead Depth\tAltReadDepth\tConsequence\tSift\tPolyPhen\tHGVSc\tHGVSp\tdbSNPID\tAlleleFreqGlobalMinor\tGlobalMinorAllele\tAlleleFreqAmr\tAlleleFreqAsn\tAlleleFreqAf\tAlleleFreqEur\tClinVarAccession\tClinVarAlleleType\tClinVarSignificance\tcosmicV72' /home/environments/$environmentID/protonAnalysis/$runName/variantCaller_out."$callerID"/$sampleName/TSVC_variants.final.txt

echo "processing amplicon coverage file"
ampliconFile=$(ls $ampliconFolder/$sampleName/*.amplicon.cov.xls)
ln -s $ampliconFile /home/environments/$environmentID/protonAnalysis/$runName/coverageAnalysis_out."$coverageID"/$sampleName/
awk -F'\t' -v OFS='\t' '{if($10 < 100) print $4,$10}' $ampliconFile > /home/environments/$environmentID/protonAnalysis/$runName/coverageAnalysis_out."$coverageID"/$sampleName/amplicon.lessThan100.txt
done

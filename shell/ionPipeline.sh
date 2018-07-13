if [ $# -eq 0 ]
then
	echo "Usage: ionPipeline.sh"
	echo "-r runID: three digits, last number of the run"
	echo "-s sampleID"
	echo "-c coverageID: coverage analysis"
	echo "-v callerID: variant caller"
	echo "-i instrumentID"
	echo "-e excluded amplicon list, one amplicon name per line (optional)"
	echo "-a amplicon bed file (optional, required if -e is specified. Only variants within ranges specified in this file will be analysed) "
	echo "-n environmentID"
	echo "-q queueID"
	exit
fi

if test $# -gt 0
	then
	while getopts :r:s:c:v:i:e:a:n:q: opt
	do
	case $opt in
		r)
			runID=$OPTARG
			;;
		s)
		  sampleID=$OPTARG
			;;
		c)
			coverageID=$OPTARG
			;;
		v)
			callerID=$OPTARG
			;;
		i)
			instrumentID=$OPTARG
			;;
		e)
			excluded=$OPTARG
			;;
		a)
			ampliconRef=$OPTARG
			;;
		n)
		  environmentID=$OPTARG
		  ;;
		q)
			  queueID=$OPTARG
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


echo "running pipeline for queue: $queueID"

###check input parameter for correctness
if [ $instrumentID != "pgm" ] && [ $instrumentID != "proton" ]
then
	echo "Error: Only pgm or proton are valid input for -i"
exit
fi

if [ -z $ampliconRef ] && [ ! -z $excluded ]
then
	echo "Error: amplicon file must be specified if the excluded list is specified"
exit
fi


if [ ! -f $ampliconRef ]
then
	echo "Error: reference amplicon file not found"
exit
fi

if [ ! -z $excluded ] && [ ! -f $excluded ]
then
	echo "Error: excluded amplicon list not found"
exit
fi

if [ -z $runID ] || [ -z $sampleID ] || [ -z $coverageID ] || [ -z $callerID ] || [ -z $instrumentID ] || [ -z $environmentID ]
then
	echo "Error: Please input required parameters-"
	echo "-r runID: three digits, last number of the run"
	echo "-s sampleID"
	echo "-c coverageID: coverage analysis"
	echo "-v callerID: variant caller"
	echo "-i instrumentID"
	echo "-e environmentID"
	exit
fi


variantFolder=$(ls -d /home/$instrumentID/*$runID/plugin_out/variantCaller_out."$callerID")
ampliconFolder=$(ls -d /home/$instrumentID/*"$runID"/plugin_out/coverageAnalysis_out."$coverageID")
runFolder=$(ls -d /home/$instrumentID/*$runID)
runName=${runFolder##*/}

#make output directory

if [ ! -d /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/variantCaller_out."$callerID" ]
then
	mkdir /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/variantCaller_out."$callerID"
fi

chmod 777 /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/variantCaller_out."$callerID"

if [ ! -d /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/coverageAnalysis_out."$coverageID" ]
then
	mkdir /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/coverageAnalysis_out."$coverageID"
fi

chmod 777 /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/coverageAnalysis_out."$coverageID"

echo "Runfolder is $runFolder"

##get runDate information##
declare -A months
months=( ["Jan"]="01" ["Feb"]="02" ["Mar"]="03" ["Apr"]="04" ["May"]="05" ["Jun"]="06" ["Jul"]="07" ["Aug"]="08" ["Sep"]="09" ["Oct"]="10" ["Nov"]="11" ["Dec"]="12" )
if [ -f /home/"$instrumentID"/*"$runID"/InitLog.txt ]
then
	runDate=$(head -n 1 /home/"$instrumentID"/*"$runID"/InitLog.txt)
	year1=$(echo $runDate |cut -d ' ' -f 5)
	year=${year1%:}
	day=$(echo $runDate |cut -d ' ' -f 3)
	monthWord=$(echo $runDate |cut -d ' ' -f 2)
	month=${months["$monthWord"]}
	date=$year-$month-$day
	echo $date > /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/runDate.txt
else
	echo "Warning: InitLog.txt file not found, run date will not be entered"
fi

# for file in $variantFolder/IonXpress_*/TSVC_variants.vcf
# for file in $variantFolder/$sampleID/TSVC_variants.vcf
# do
## note: removed for loop, passed parameter sampleID is same as sampleName which
## is previously used for running all the downstream pipelines
## do not confuse with similar variable sample_ID which is only used once for horizon DNA

file=$variantFolder/$sampleID/TSVC_variants.vcf
echo "Processing $file"
sampleFolder=${file%/*}
sampleName=${sampleFolder##*/}
if [ ! -d /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/variantCaller_out."$callerID" ]
then
	mkdir /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/variantCaller_out."$callerID"
fi
if [ ! -d /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/coverageAnalysis_out."$coverageID" ]
then
	mkdir /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/coverageAnalysis_out."$coverageID"
fi

chmod 777 /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/variantCaller_out."$callerID"
chmod 777 /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/coverageAnalysis_out."$coverageID"

echo "filter against amplicon"
echo " sampleID is $sampleID"
sample_ID=$(grep "^#CHROM" $file |cut -f 10)
echo " sample Folder is $sampleFolder"
echo " sampleName is $sampleName"
echo " sample_ID is $sample_ID"


echo "----------> File is $file"
echo "----------> File is $ampliconRef"

bash /home/pipelines/master/shell/updatepipelineStatus.sh -q $queueID -s bedtools

if [ ! -z $ampliconRef ]
then
	/opt/software/bedtools-2.17.0/bin/bedtools intersect -u -a $file -b $ampliconRef > /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/variantCaller_out."$callerID"/TSVC_variants.filter.vcf
else
	ln -s $file /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/variantCaller_out."$callerID"/TSVC_variants.filter.vcf
fi


bash /home/pipelines/master/shell/updatepipelineStatus.sh -q $queueID -s splitVcf

echo "split multiple alt alleles into different lines"
python /home/pipelines/master/python/splitVcf.py \
-I /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/variantCaller_out."$callerID"/TSVC_variants.filter.vcf \
-f FAO,FDP,AF \
-o /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/variantCaller_out."$callerID"/TSVC_variants.split.vcf

bash /home/pipelines/master/shell/updatepipelineStatus.sh -q $queueID -s vep

echo "running VEP"
/home/pipelines/master/perl/bin/perl \
/opt/vep/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \
-i /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/variantCaller_out."$callerID"/TSVC_variants.split.vcf \
-o /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/variantCaller_out."$callerID"/TSVC_variants.split.vep.vcf \
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
parseIonNewVarView \
-I /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/variantCaller_out."$callerID"/TSVC_variants.split.vep.vcf \
-o /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/variantCaller_out."$callerID"/TSVC_variants.split.vep.parse.newVarView.txt

echo "filter VEP results"
shopt -s nocasematch
if [[ $sample_ID =~ horizon ]]
then
	awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 1 && $10 != "null" && $11 >=100 && $11 != "null") print}' /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/variantCaller_out."$callerID"/TSVC_variants.split.vep.parse.newVarView.txt \
	> /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/variantCaller_out."$callerID"/TSVC_variants.split.vep.parse.newVarView.filter.txt
else
	awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 10 && $10 != "null" && $11 >=100 && $11 != "null") print}' /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/variantCaller_out."$callerID"/TSVC_variants.split.vep.parse.newVarView.txt \
	> /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/variantCaller_out."$callerID"/TSVC_variants.split.vep.parse.newVarView.filter.txt
fi


sed -i '1iGene\texon\tchr\tpos\tref\talt\tClassification\tType\tQuality\tAltVariantFreq\tRead Depth\tAltReadDepth\tConsequence\tSift\tPolyPhen\tHGVSc\tHGVSp\tdbSNPID\tpubmed' /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/variantCaller_out."$callerID"/TSVC_variants.split.vep.parse.newVarView.filter.txt

echo "processing amplicon coverage file"
ampliconFile=$(ls $ampliconFolder/$sampleName/*.amplicon.cov.xls)
if [ ! -z $excluded ]
then
	grep -v -f $excluded $ampliconFile > /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/coverageAnalysis_out."$coverageID"/amplicon.filter.txt
else
	ln -s "$ampliconFile" /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/coverageAnalysis_out."$coverageID"/amplicon.filter.txt
fi

awk -F'\t' -v OFS='\t' '{if($10 < 100) print $4,$10}' /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/coverageAnalysis_out."$coverageID"/amplicon.filter.txt \
> /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/coverageAnalysis_out."$coverageID"/amplicon.lessThan100.txt


bash /home/pipelines/master/shell/updatepipelineStatus.sh -q $queueID -s runCompleted

echo "Pipeline Finished!"

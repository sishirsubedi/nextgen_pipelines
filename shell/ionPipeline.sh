if [ $# -eq 0 ]
then
	echo "Usage: ionPipeline.sh"
	echo "-r runID: three digits, last number of the run"
	echo "-s sampleName"
	echo "-c coverageID: coverage analysis"
	echo "-v callerID: variant caller"
	echo "-i instrument"
	echo "-e excluded amplicon list, one amplicon name per line (optional)"
	echo "-a amplicon bed file (optional, required if -e is specified. Only variants within ranges specified in this file will be analysed) "
	echo "-n environment"
	echo "-q queueID"
	echo "-u user"
	echo "-p password"
	exit
fi


if test $# -gt 0
	then
	while getopts :r:s:c:v:i:e:a:n:q:u:p: opt
	do
	case $opt in
		r)
			runID=$OPTARG
			;;
		s)
		  sampleName=$OPTARG
			;;
		c)
			coverageID=$OPTARG
			;;
		v)
			callerID=$OPTARG
			;;
		i)
			instrument=$OPTARG
			;;
		e)
			excluded=$OPTARG
			;;
		a)
			ampliconRef=$OPTARG
			;;
		n)
		  environment=$OPTARG
		  ;;
		q)
			queueID=$OPTARG
			  ;;
		u)
			user=$OPTARG
				;;
		p)
			 password=$OPTARG
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



function updateStatus() {
user=$4
password=$5
database=$3
insertstatement="INSERT INTO pipelineStatus (queueID, plStatus, timeUpdated) VALUES ('$1','$2',now());"
mysql --user="$user" --password="$password" --database="$database" --execute="$insertstatement"
}

echo " $currentdate    INFO  -  running ionPipeline for queue: $queueID"

###check input parameter for correctness
if [ $instrument != "pgm" ] && [ $instrument != "proton" ]
then
	echo "Error: Only pgm or proton are valid input for -i"
	updateStatus "$queueID" "ERROR:instrument" "$environment" "$user"  "$password"
  exit
fi

if [ -z $ampliconRef ] && [ ! -z $excluded ]
then
	echo "Error: amplicon file must be specified if the excluded list is specified"
  updateStatus "$queueID" "ERROR:amp_exc" "$environment" "$user"  "$password"
	exit
fi


if [ ! -f $ampliconRef ]
then
	echo "Error: reference amplicon file not found"
	updateStatus "$queueID" "ERROR:amplicon" "$environment" "$user"  "$password"
  exit
fi

if [ ! -z $excluded ] && [ ! -f $excluded ]
then
	echo "Error: excluded amplicon list not found"
  updateStatus "$queueID" "ERROR:excluded" "$environment" "$user"  "$password"
  exit
fi

if [ -z $runID ] || [ -z $sampleName ] || [ -z $coverageID ] || [ -z $callerID ] || [ -z $instrument ] || [ -z $environment ]
then
	echo "Error: Please input required parameters-"
	echo "-r runID: three digits, last number of the run"
	echo "-s sampleName"
	echo "-c coverageID: coverage analysis"
	echo "-v callerID: variant caller"
	echo "-i instrument"
	echo "-e environment"
	updateStatus "$queueID" "ERROR:parameters" "$environment" "$user"  "$password"
	exit
fi


variantFolder=$(ls -d /home/$instrument/*$runID/plugin_out/"$callerID")
ampliconFolder=$(ls -d /home/$instrument/*"$runID"/plugin_out/"$coverageID")
runFolder=$(ls -d /home/$instrument/*$runID)
runName=${runFolder##*/}

#make output directory

if [ ! -d /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID" ]
then
	mkdir /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"
fi

chmod 777 /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"

if [ ! -d /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$coverageID" ]
then
	mkdir /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$coverageID"
fi

chmod 777 /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$coverageID"

echo " $currentdate    INFO  -  Runfolder is $runFolder"

##get runDate information##
declare -A months
months=( ["Jan"]="01" ["Feb"]="02" ["Mar"]="03" ["Apr"]="04" ["May"]="05" ["Jun"]="06" ["Jul"]="07" ["Aug"]="08" ["Sep"]="09" ["Oct"]="10" ["Nov"]="11" ["Dec"]="12" )
if [ -f /home/"$instrument"/*"$runID"/InitLog.txt ]
then
	runDate=$(head -n 1 /home/"$instrument"/*"$runID"/InitLog.txt)
	year1=$(echo $runDate |cut -d ' ' -f 5)
	year=${year1%:}
	day=$(echo $runDate |cut -d ' ' -f 3)
	monthWord=$(echo $runDate |cut -d ' ' -f 2)
	month=${months["$monthWord"]}
	date=$year-$month-$day
	echo $date > /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/runDate.txt
else
	echo "Warning: InitLog.txt file not found, run date will not be entered"
fi

# for file in $variantFolder/IonXpress_*/TSVC_variants.vcf
# for file in $variantFolder/$sampleName/TSVC_variants.vcf
# do
## note: removed for loop, passed parameter sampleName is same as sampleName which
## is previously used for running all the downstream pipelines
## do not confuse with similar variable sample_ID which is only used once for horizon DNA

file=$variantFolder/$sampleName/TSVC_variants.vcf
echo " $currentdate    INFO  -  Processing $file"
sampleFolder=${file%/*}
sampleName=${sampleFolder##*/}
if [ ! -d /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID" ]
then
	mkdir /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"
fi
if [ ! -d /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$coverageID" ]
then
	mkdir /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$coverageID"
fi

chmod 777 /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"
chmod 777 /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$coverageID"

sample_ID=$(grep "^#CHROM" $file |cut -f 10)
updateStatus "$queueID" "bedtools" "$environment" "$user"  "$password"

if [ ! -z $ampliconRef ]
then
	/opt/software/bedtools-2.17.0/bin/bedtools intersect -u -a $file -b $ampliconRef > /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"/TSVC_variants.filter.vcf
else
	ln -s $file /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"/TSVC_variants.filter.vcf
fi


if [ ! -f /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"/TSVC_variants.filter.vcf ]
then
	echo "Error: bedtools"
	updateStatus "$queueID" "ERROR:bedtools" "$environment" "$user"  "$password"
exit
fi

updateStatus "$queueID" "splitVcf" "$environment" "$user"  "$password"

echo " $currentdate    INFO  -  split multiple alt alleles into different lines"
python /var/pipelines_"$environment"/python/splitVcf.py \
-I /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"/TSVC_variants.filter.vcf \
-f FAO,FDP,AF \
-o /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"/TSVC_variants.split.vcf


if [ ! -f /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"/TSVC_variants.split.vcf ]
then
	echo "Error: splitVcf"
	updateStatus "$queueID" "ERROR:splitVcf" "$environment" "$user"  "$password"
  exit
fi


updateStatus "$queueID" "VEP" "$environment" "$user"  "$password"

echo " $currentdate    INFO  -  running VEP"
/opt/perl/bin/perl \
/opt/vep/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \
-i /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"/TSVC_variants.split.vcf \
-o /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"/TSVC_variants.split.vep.vcf \
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



echo " $currentdate    INFO  -  parse VEP results"
python /var/pipelines_"$environment"/python/parseVEP.py \
parseIonNewVarView \
-I /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"/TSVC_variants.split.vep.vcf \
-o /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"/TSVC_variants.split.vep.parse.newVarView.txt

if [ ! -f /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"/TSVC_variants.split.vep.parse.newVarView.txt ]
then
	echo "Error: vep"
	updateStatus "$queueID" "ERROR:VEP" "$environment" "$user"  "$password"
  exit
fi



echo " $currentdate    INFO  -  filter VEP results"
shopt -s nocasematch
if [[ $sample_ID =~ horizon ]]
then
	awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 1 && $10 != "null" && $11 >=100 && $11 != "null") print}' /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"/TSVC_variants.split.vep.parse.newVarView.txt \
	> /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"/TSVC_variants.split.vep.parse.newVarView.filter.txt
else
	awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 10 && $10 != "null" && $11 >=100 && $11 != "null") print}' /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"/TSVC_variants.split.vep.parse.newVarView.txt \
	> /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"/TSVC_variants.split.vep.parse.newVarView.filter.txt
fi


echo " $currentdate    INFO  -  running sed"
sed -i '1iGene\texon\tchr\tpos\tref\talt\tClassification\tType\tQuality\tAltVariantFreq\tRead Depth\tAltReadDepth\tConsequence\tSift\tPolyPhen\tHGVSc\tHGVSp\tdbSNPID\tpubmed' /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"/TSVC_variants.split.vep.parse.newVarView.txt
sed -i '1iGene\texon\tchr\tpos\tref\talt\tClassification\tType\tQuality\tAltVariantFreq\tRead Depth\tAltReadDepth\tConsequence\tSift\tPolyPhen\tHGVSc\tHGVSp\tdbSNPID\tpubmed' /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"/TSVC_variants.split.vep.parse.newVarView.filter.txt

echo " $currentdate    INFO  -  processing amplicon coverage file"
ampliconFile=$(ls $ampliconFolder/$sampleName/*.amplicon.cov.xls)
if [ ! -z $excluded ]
then
	grep -v -f $excluded $ampliconFile > /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$coverageID"/amplicon.filter.txt
else
	ln -s "$ampliconFile" /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$coverageID"/amplicon.filter.txt
fi

awk -F'\t' -v OFS='\t' '{if($10 < 100) print $4,$10}' /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$coverageID"/amplicon.filter.txt \
> /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$coverageID"/amplicon.lessThan100.txt


updateStatus "$queueID" "UpdateDatabase" "$environment" "$user"  "$password"


#### wait for database update

sleep 1s

################################################################################
# updating analysis results
################################################################################

# update analysis
analysischeck_statement="select queueID,plStatus from pipelineStatus where plStatus='UpdateDatabase' and queueID=$queueID ;"
while  read -r queueID plStatus;
do

	echo " $currentdate    INFO  -  Running analysis"

  newStatus='UpdatingDatabase'
  updateanalysis_statement="update pipelineStatus set plStatus='$newStatus' where queueID=$queueID and plStatus='$plStatus'"
  mysql --user="$user" --password="$password" --database="$environment" --execute="$updateanalysis_statement"

  bash /var/pipelines_"$environment"/shell/runAnalysis.sh -q $queueID -e $environment -u $user -p $password

done< <(mysql --user="$user" --password="$password" --database="$environment" --execute="$analysischeck_statement" -N)

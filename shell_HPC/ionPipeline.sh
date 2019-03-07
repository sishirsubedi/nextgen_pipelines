##############################################################################
#
# Houston Methodist Hospital
# Molecular Diagnostic
#
#Description:
#This script runs ion pipeline.
##############################################################################

#!/bin/bash

################################################################################
# assign Variables
################################################################################

while getopts :r:s:c:v:i:e:a:n:q:u:p: opt; do
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

################################################################################
# functions
################################################################################

log()
{
 MESSAGE=$1
 TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
 SCRIPT=$( basename $0 )
 echo " [ $TIMESTAMP ] [ $SCRIPT ] : $MESSAGE "
}

function updateStatus() {
user=$4
password=$5
database="ngs_$3"
insertstatement="INSERT INTO pipelineStatus (queueID, plStatus, timeUpdated) VALUES ('$1','$2',now());"
mysql --host="hhplabngsp01" --user="$user" --password="$password" --database="$database" --execute="$insertstatement"
}

################################################################################
#
################################################################################

log " Running ionPipeline for :
environment : $environment
instrument : $instrument
assay : $assay
runID : $runID
sampleName : $sampleName
coverageID : $coverageID
callerID : $callerID
queueID : $queueID"


if [ $instrument != "proton" ]
then
	log "Error: Only pgm or proton are valid input for -i"
	updateStatus "$queueID" "ERROR:instrument" "$environment" "$user"  "$password"
  exit
fi

if [ -z $ampliconRef ] && [ ! -z $excluded ]
then
	log "Error: amplicon file must be specified if the excluded list is specified"
  updateStatus "$queueID" "ERROR:amp_exc" "$environment" "$user"  "$password"
	exit
fi


if [ ! -f $ampliconRef ]
then
	log "Error: reference amplicon file not found"
	updateStatus "$queueID" "ERROR:amplicon" "$environment" "$user"  "$password"
  exit
fi

if [ ! -z $excluded ] && [ ! -f $excluded ]
then
	log "Error: excluded amplicon list not found"
  updateStatus "$queueID" "ERROR:excluded" "$environment" "$user"  "$password"
  exit
fi

if [ -z $runID ] || [ -z $sampleName ] || [ -z $coverageID ] || [ -z $callerID ] || [ -z $instrument ] || [ -z $environment ]
then
	log "Error: Please input required parameters-
	-r runID: three digits, last number of the run
	-s sampleName
	-c coverageID: coverage analysis
	-v callerID: variant caller
	-i instrument
	-e environment"
	updateStatus "$queueID" "ERROR:parameters" "$environment" "$user"  "$password"
	exit
fi

################################################################################
# initialize variables
################################################################################

variantFolder=$(ls -d /home/${instrument}/*${runID}/plugin_out/${callerID})
ampliconFolder=$(ls -d /home/${instrument}/*${runID}/plugin_out/${coverageID})
runFolder=$(ls -d /home/$instrument/*$runID)
runName=${runFolder##*/}
ENV_HOME="/home/environments/ngs_${environment}/${instrument}Analysis/$runName/$sampleName/"
SCRIPT_HOME="/home/pipelines/ngs_${environment}/"
SAMPLE_FILE="$variantFolder/${sampleName}/TSVC_variants.vcf"

log " $currentdate    INFO  -  Processing $SAMPLE_FILE"

################################################################################
#
################################################################################

# if [ ! -d ${ENV_HOME}$callerID ]
# then
# 	mkdir ${ENV_HOME}$callerID
# fi
# chmod 775 ${ENV_HOME}$callerID
#
# if [ ! -d ${ENV_HOME}$coverageID ]
# then
# 	mkdir ${ENV_HOME}$coverageID
# fi
# chmod 775 ${ENV_HOME}$coverageID

##get runDate information
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
	echo $date > ${ENV_HOME}runDate.txt
else
	log "Warning: InitLog.txt SAMPLE_FILE not found, run date will not be entered"
fi

################################################################################
# bedtools
################################################################################

updateStatus "$queueID" "bedtools" "$environment" "$user"  "$password"

if [ ! -z $ampliconRef ]
then
	/opt/software/bedtools-2.17.0/bin/bedtools intersect -u -a $SAMPLE_FILE -b $ampliconRef > ${ENV_HOME}${callerID}/TSVC_variants.filter.vcf
else
	ln -s $SAMPLE_FILE ${ENV_HOME}${callerID}/TSVC_variants.filter.vcf
fi


if [ ! -f ${ENV_HOME}${callerID}/TSVC_variants.filter.vcf ]
then
	log "Error: bedtools"
	updateStatus "$queueID" "ERROR:FILE-NOT-FOUND:TSVC_variants" "$environment" "$user"  "$password"
exit
fi

updateStatus "$queueID" "splitVcf" "$environment" "$user"  "$password"

log " $currentdate    INFO  -  split multiple alt alleles into different lines"
python ${SCRIPT_HOME}python/splitVcf.py \
-I ${ENV_HOME}${callerID}/TSVC_variants.filter.vcf \
-f FAO,FDP,AF \
-o ${ENV_HOME}${callerID}/TSVC_variants.filter.split.vcf


if [ ! -f ${ENV_HOME}${callerID}/TSVC_variants.filter.split.vcf ]
then
	log "Error: splitVcf"
	updateStatus "$queueID" "ERROR:FILE-NOT-FOUND:splitVcf" "$environment" "$user"  "$password"
  exit
fi


updateStatus "$queueID" "VEP" "$environment" "$user"  "$password"

log " $currentdate    INFO  -  running VEP"
/opt/perl/bin/perl /opt/vep/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \
-i ${ENV_HOME}${callerID}/TSVC_variants.filter.split.vcf \
-o ${ENV_HOME}${callerID}/TSVC_variants.filter.split.vep.vcf \
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


log " $currentdate    INFO  -  parse VEP results"
python ${SCRIPT_HOME}/python/parseVEP.py \
parseIonNewVarView \
-I ${ENV_HOME}${callerID}/TSVC_variants.filter.split.vep.vcf \
-o ${ENV_HOME}${callerID}/TSVC_variants.filter.split.vep.parse.txt

if [ ! -f ${ENV_HOME}${callerID}/TSVC_variants.filter.split.vep.parse.txt ]
then
	log "Error: vep"
	updateStatus "$queueID" "ERROR:VEP" "$environment" "$user"  "$password"
  exit
fi


log " $currentdate    INFO  -  filter VEP results"
shopt -s nocasematch
if [[ $sampleName =~ horizon ]]
then
	awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 1 && $10 != "null" && $11 >=100 && $11 != "null") print}' ${ENV_HOME}${callerID}/TSVC_variants.filter.split.vep.parse.txt \
	> ${ENV_HOME}${callerID}/TSVC_variants.filter.split.vep.parse.filter2.txt
else
	awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 10 && $10 != "null" && $11 >=100 && $11 != "null") print}' ${ENV_HOME}${callerID}/TSVC_variants.filter.split.vep.parse.txt \
	> ${ENV_HOME}${callerID}/TSVC_variants.filter.split.vep.parse.filter2.txt
fi

# log " $currentdate    INFO  -  running sed"
# sed -i '1iGene\texon\tchr\tpos\tref\talt\tClassification\tType\tQuality\tAltVariantFreq\tRead Depth\tAltReadDepth\tConsequence\tSift\tPolyPhen\tHGVSc\tHGVSp\tdbSNPID\tpubmed' /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"/TSVC_variants.split.vep.parse.newVarView.txt
# sed -i '1iGene\texon\tchr\tpos\tref\talt\tClassification\tType\tQuality\tAltVariantFreq\tRead Depth\tAltReadDepth\tConsequence\tSift\tPolyPhen\tHGVSc\tHGVSp\tdbSNPID\tpubmed' /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$callerID"/TSVC_variants.split.vep.parse.newVarView.filter.txt

log " $currentdate    INFO  -  processing amplicon coverage file"
ampliconFile=$(ls $ampliconFolder/${sampleName}/*.amplicon.cov.xls)
if [ ! -z $excluded ]
then
	grep -v -f $excluded $ampliconFile > ${ENV_HOME}${coverageID}/amplicon.filter.txt
else
	ln -s "$ampliconFile" ${ENV_HOME}${coverageID}/amplicon.filter.txt
fi

awk -F'\t' -v OFS='\t' '{if($10 < 100) print $4,$10}' ${ENV_HOME}${coverageID}/amplicon.filter.txt \
> ${ENV_HOME}${coverageID}/amplicon.lessThan100.txt


updateStatus "$queueID" "UpdateDatabase" "$environment" "$user"  "$password"


################################################################################
# updating analysis results
################################################################################

# update analysis
# analysischeck_statement="select queueID,plStatus from pipelineStatus where plStatus='UpdateDatabase' and queueID=$queueID ;"
# while  read -r queueID plStatus;
# do
#
# 	log " $currentdate    INFO  -  Running analysis"
#
#   newStatus='UpdatingDatabase'
#   updateanalysis_statement="update pipelineStatus set plStatus='$newStatus' where queueID=$queueID and plStatus='$plStatus'"
#   mysql --host="hhplabngsp01" --user="$user" --password="$password" --database="$environment" --execute="$updateanalysis_statement"
#
#   bash /var/pipelines_"$environment"/shell/runAnalysis.sh -q $queueID -e $environment -u $user -p $password
#
# done< <(mysql --host="hhplabngsp01" --user="$user" --password="$password" --database="$environment" --execute="$analysischeck_statement" -N)

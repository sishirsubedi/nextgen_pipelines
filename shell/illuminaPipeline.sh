#!/bin/bash

if [ $# -eq 0 ]
then
echo "Usage: illuminaPipeline.sh"
echo "-r runID: three digits, last number of the run"
echo "-s sampleName"
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
	while getopts :r:s:i:e:a:n:q:u:p: opt
	do
	case $opt in
	r)
		runID=$OPTARG
		;;
	s)
	  sampleName=$OPTARG
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


if [ $instrument == "miseq" ]
then
		if [ ! -d /home/$instrument/*_"$runID"_*/Data/Intensities/BaseCalls/Alignment/ ]
		then
		echo "Error: result folder for run ID:$runID not found"
		updateStatus "$queueID" "ERROR:alignmentFolder" "$environment" "$user"  "$password"

		exit
		fi

		if [ ! -f /home/$instrument/*_"$runID"_*/Data/Intensities/BaseCalls/Alignment/AmpliconCoverage_M1.tsv ]
		then
		echo "Error: amplicon file for run ID:$runID not found"
		updateStatus "$queueID" "ERROR:ampliconcoverage_file" "$environment" "$user"  "$password"

		exit
	  fi
elif [ $instrument == "nextseq" ]
then
		if [ ! -d /home/$instrument/*_"$runID"_*/ ]
		then
		echo "Error: result folder for run ID:$runID not found"
		updateStatus "$queueID" "ERROR:runID_file" "$environment" "$user"  "$password"
		exit
	  fi
fi

if [ -z $instrument ]
then
echo "Error: instrument must be specificed with -i option"
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
updateStatus "$queueID" "ERROR:amplicon"  "$environment" "$user"  "$password"

exit
fi

if [ ! -z $excluded ] && [ ! -f $excluded ]
then
echo "Error: excluded amplicon list not found"
updateStatus "$queueID" "ERROR:excluded"  "$environment" "$user"  "$password"

exit
fi

if [ -z $runID ] || [ -z $instrument ]
then
echo "Error: -r runID and -i instrument are the required parameters"
updateStatus "$queueID" "ERROR:run_ins"  "$environment" "$user"  "$password"

exit
fi

if [ $instrument == "miseq" ]
then

	ampliconFile=$(ls /home/$instrument/*_"$runID"_*/Data/Intensities/BaseCalls/Alignment/AmpliconCoverage_M1.tsv)
	runFolder=$(ls -d /home/$instrument/*_"$runID"_*/Data/Intensities/BaseCalls/Alignment/)
	post=${runFolder##/home/$instrument/}
	runName=${post%%/Data/Intensities/BaseCalls/Alignment*}


	##get runDate information##
	dateString=${runName%%_*}
	year=20${dateString:0:2}
	month=${dateString:2:2}
	day=${dateString:4:2}
	dateString=$year-$month-$day

	echo $dateString > /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/runDate.txt


	###filter amplicon file
	if [ ! -z $excluded ]
	then
		grep -v -f $excluded $ampliconFile > /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.filtered.txt
	else
		ln -s $ampliconFile /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.filtered.txt
	fi

	###Start processing file
	for file in $runFolder/*.vcf
	do
		if  [[ ${file} =~ $sampleName ]]
	  then
			break
		fi
	done
	#
  echo "working file is - " $file


	# if [[ $file =~ .*_S[0-9]+\.vcf ]]
	# #if the vcf is in the raw machine output form
	# then
	#
	# fileName=${file##*/}
	# sampleName=${fileName%%_*}
	# fi

  updateStatus "$queueID" "bedtools"  "$environment"  "$user"  "$password"
  echo   "ampliconref is "$ampliconRef
	##filter vcf against ampliconRef
	if [ ! -z $ampliconRef ]   #string is null, that is, has zero length
	then
	/opt/software/bedtools-2.17.0/bin/bedtools intersect -u -a $file -b $ampliconRef > /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.vcf
	else
	ln -s $file /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.vcf
	fi

	if [ ! -f /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.vcf ]
	then
		echo "Error:bedtools"
		updateStatus "$queueID" "ERROR:bedtools"  "$environment" "$user"  "$password"
	exit
	fi

	updateStatus "$queueID" "VEP"  "$environment" "$user"  "$password"

	echo " ''$currentdate''    INFO  - running VEP"
	/opt/perl/bin/perl \
	/opt/vep/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \
	-i /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.vcf \
	-o /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.vep.vcf \
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

	if [ ! -f /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.vep.vcf ]
	then
		echo "Error:vep"
		updateStatus "$queueID" "ERROR:VEP"  "$environment" "$user"  "$password"
	exit
	fi


  updateStatus "$queueID" "parseVEP"  "$environment" "$user"  "$password"

	echo " ''$currentdate''    INFO  - parse VEP results for new varView"
	python /var/pipelines_"$environment"/python/parseVEP.py \
	parseIllumina \
	-I /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.vep.vcf \
	-o /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.vep.parse.txt

	echo " ''$currentdate''    INFO  - filter VEP results"
	shopt -s nocasematch
	if [[ $sampleName =~ horizon ]]
	then
	awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 1 && $11 >=100) print}' /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.vep.parse.txt \
	> /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.vep.parse.filter.txt
	else
	awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 10 && $11 >=100) print}' /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.vep.parse.txt \
	> /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.vep.parse.filter.txt
	fi


	sed -i '1iGene\texon\tchr\tpos\tref\talt\tClassification\tType\tQuality\tAltVariantFreq\tRead Depth\tAltReadDepth\tConsequence\tSift\tPolyPhen\tHGVSc\tHGVSp\tdbSNPID\tpubmed' /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.vep.parse.filter.txt

	##split amplicon file

	sampleCol=$(head -n 1 /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.filtered.txt | awk -F"\t" -v sample=$sampleName '{for (i=1; i<= NF; i++) if($i==sample) print i}')

	if [ $sampleCol = "" ]
		then
			echo "amplicon info not found in amplicon file for sample:$sampleName"
	fi

	cut -f 1,$sampleCol /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.filtered.txt > /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.txt
	awk '$2 < 100' /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.txt > /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.lessThan100.txt


  updateStatus  "$queueID" "UpdateDatabase"  "$environment" "$user"  "$password"

elif  [ $instrument == "nextseq" ]
then

	runFolder=$(ls -d /home/$instrument/*_"$runID"_*)
	post=${runFolder##/home/$instrument/}
	runName=${post%%/Data/Intensities/BaseCalls/Alignment*}


	if [ ! -d /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis ]
	then
		mkdir /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis
	fi
	chmod 777 /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis

	##get runDate information##
	dateString=${runName%%_*}
	year=20${dateString:0:2}
	month=${dateString:2:2}
	day=${dateString:4:2}
	dateString=$year-$month-$day

	echo $dateString > /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/runDate.txt


	###Start processing file
	#for file in $runFolder/*.vcf
	for file in /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantCaller/"$sampleName".comb.vcf
	do
		echo $file
		if  [[ ${file} =~ $sampleName ]]
		then
			break
		fi
	done

	echo "file is $file"

	if [[ $file =~ .comb.vcf ]]
		##if the vcf is in the raw machine output form
	then
		echo "Processing $file "
		fileName=${file##*/}
		#sampleName=${fileName%%_*}
		sampleName=${fileName%%.comb*}

		echo "fileName is $fileName"
		echo "sampleName is $sampleName"


    updateStatus "$queueID" "bedtools"  "$environment"  "$user"  "$password"

			##filter vcf against ampliconRef
		if [ ! -z $ampliconRef ]
			then
			/opt/software/bedtools-2.17.0/bin/bedtools intersect -u -a $file -b $ampliconRef > /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.vcf
		else
			ln -s $file /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.vcf         #create symbolic link
		fi

		if [ ! -f /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.vcf ]
	 	then
	 		echo "Error:bedtools"
	 		updateStatus "$queueID" "ERROR:bedtools"  "$environment" "$user"  "$password"
	 	exit
	 	fi



    updateStatus "$queueID" "VEP"   "$environment"  "$user"  "$password"

		echo " ''$currentdate''    INFO  - running VEP"
		/opt/perl/bin/perl \
		/opt/vep/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \
			-i /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.vcf \
			-o /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.vep.vcf \
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


			echo " ''$currentdate''    INFO  - parse VEP results for new varView"
			python /var/pipelines_"$environment"/python/parseVEP.py \
			parseIlluminaNextseq \
			-I /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.vep.vcf \
			-o /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.vep.parse.txt

			if [ ! -f /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.vep.parse.txt ]
			then
				echo "Error:vep"
				updateStatus "$queueID" "ERROR:VEP"  "$environment"  "$user"  "$password"
			exit
			fi



			echo " ''$currentdate''    INFO  - filter VEP results"
			shopt -s nocasematch
			if [[ $sampleName =~ horizon ]]
				then
				awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 1 && $11 >=100) print}' /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.vep.parse.txt \
				> /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.vep.parse.filter.txt
			else
				awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 10 && $11 >=100) print}' /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.vep.parse.txt \
				> /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.vep.parse.filter.txt
			fi

			echo " ''$currentdate''    INFO  - running sed"
			sed -i '1iGene\texon\tchr\tpos\tref\talt\tClassification\tType\tQuality\tAltVariantFreq\tRead Depth\tAltReadDepth\tConsequence\tSift\tPolyPhen\tHGVSc\tHGVSp\tdbSNPID\tpubmed' /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.vep.parse.txt
      sed -i '1iGene\texon\tchr\tpos\tref\talt\tClassification\tType\tQuality\tAltVariantFreq\tRead Depth\tAltReadDepth\tConsequence\tSift\tPolyPhen\tHGVSc\tHGVSp\tdbSNPID\tpubmed' /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.vep.parse.filter.txt

      updateStatus "$queueID" "samtools"  "$environment"  "$user"  "$password"

			echo " ''$currentdate''    INFO  - generating "$sampleName" ampliconFile"
			/opt/samtools-1.4/samtools-1.4/samtools bedcov  /doc/ref/Heme/trusight-myeloid-amplicon-track.excludedNew.bed /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantCaller/$sampleName.sort.bam | awk ' {print $4,"\t",int($13/($8-$7))} ' > /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.samtools.coverageDepth
			ampliconFile=$(ls /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.samtools.coverageDepth)

			###filter amplicon file
			if [ ! -z $excluded ]
			then
			grep -v -f $excluded $ampliconFile > /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.filtered.txt
			else
			ln -s $ampliconFile /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.filtered.txt
			fi



			##split amplicon file
			sampleCol=$(head -n 1 /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.filtered.txt | awk -F"\t" '{print NF; exit}')
			if [ $sampleCol != "2" ]
				then
				echo "amplicon info not found in amplicon file for sample:$sampleName"
			fi


			cut -f 1,$sampleCol /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.filtered.txt > /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.txt
			awk '$2 < 100' /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.txt > /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.lessThan100.txt

			updateStatus "$queueID" "UpdateDatabase"  "$environment"  "$user"  "$password"

	fi


fi

#### wait for database update

sleep 1

################################################################################
# updating analysis results
################################################################################

analysischeck_statement="select queueID,plStatus from pipelineStatus where plStatus='UpdateDatabase' and queueID=$queueID ;"
while  read -r queueID plStatus;
do
  newStatus='UpdatingDatabase'
  echo " ''$currentdate''    INFO  - Running analysis"
  updateanalysis_statement="update pipelineStatus set plStatus='$newStatus' where queueID=$queueID and plStatus='$plStatus'"
  mysql --user="$user" --password="$password" --database="$environment" --execute="$updateanalysis_statement"
  bash /var/pipelines_"$environment"/shell/runAnalysis.sh -q $queueID -e $environment -u $user -p $password
done< <(mysql --user="$user" --password="$password" --database="$environment" --execute="$analysischeck_statement" -N)

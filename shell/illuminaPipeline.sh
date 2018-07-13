if [ $# -eq 0 ]
then
echo "Usage: illuminaPipeline.sh"
echo "-r runID: three digits, last number of the run"
echo "-s sampleID"
echo "-i instrumentID"
echo "-e excluded amplicon list, one amplicon name per line (optional)"
echo "-a amplicon bed file (optional, required if -e is specified. Only variants within ranges specified in this file will be analysed) "
echo "-n environmentID"
echo "-q queueID"
exit
fi

if test $# -gt 0
	then
	while getopts :r:s:i:e:a:n:q: opt
	do
	case $opt in
	r)
		runID=$OPTARG
		;;
	s)
	  sampleID=$OPTARG
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


if [ $instrumentID == "miseq" ]
then
		if [ ! -d /home/$instrumentID/*_"$runID"_*/Data/Intensities/BaseCalls/Alignment/ ]
		then
		echo "Error: result folder for run ID:$runID not found"
		exit
		fi

		if [ ! -f /home/$instrumentID/*_"$runID"_*/Data/Intensities/BaseCalls/Alignment/AmpliconCoverage_M1.tsv ]
		then
		echo "Error: amplicon file for run ID:$runID not found"
		exit
	  fi
elif [ $instrumentID == "nextseq" ]
then
		if [ ! -d /home/$instrumentID/*_"$runID"_*/]
		then
		echo "Error: result folder for run ID:$runID not found"
		exit
	  fi
fi

if [ -z $instrumentID ]
then
echo "Error: instrument must be specificed with -i option"
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

if [ -z $runID ] || [ -z $instrumentID ]
then
echo "Error: -r runID and -i instrumentID are the required parameters"
exit
fi

if [ $instrumentID == "miseq" ]
then

	ampliconFile=$(ls /home/$instrumentID/*_"$runID"_*/Data/Intensities/BaseCalls/Alignment/AmpliconCoverage_M1.tsv)
	runFolder=$(ls -d /home/$instrumentID/*_"$runID"_*/Data/Intensities/BaseCalls/Alignment/)
	post=${runFolder##/home/$instrumentID/}
	runName=${post%%/Data/Intensities/BaseCalls/Alignment*}


	echo "Runfolder is $runFolder"

	##get runDate information##
	dateString=${runName%%_*}
	year=20${dateString:0:2}
	month=${dateString:2:2}
	day=${dateString:4:2}
	dateString=$year-$month-$day

	echo $dateString > /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/runDate.txt


	###filter amplicon file
	if [ ! -z $excluded ]
	then
		grep -v -f $excluded $ampliconFile > /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/amplicon.filtered.txt
	else
		ln -s $ampliconFile /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/amplicon.filtered.txt
	fi

	echo "runfolder is $runFolder"
	echo "sample id is $sampleID"
	###Start processing file
	for file in $runFolder/*.vcf
	do
		echo $file
		if  [[ ${file} =~ $sampleID ]]
	  then
			break
		fi
	done
	#
	# file=$runFolder$sampleID.vcf

	if [[ $file =~ .*_S[0-9]+\.vcf ]]
	#if the vcf is in the raw machine output form
	then
	echo "Processing $file"
	fileName=${file##*/}
	sampleName=${fileName%%_*}
	fi


	echo "----------> File is $file"
	echo "----------> File is $ampliconRef"

  bash /home/pipelines/master/shell/updatepipelineStatus.sh -q $queueID -s bedtools

	##filter vcf against ampliconRef
	if [ ! -z $ampliconRef ]
	then
	/opt/software/bedtools-2.17.0/bin/bedtools intersect -u -a $file -b $ampliconRef > /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/$sampleName.amplicon.vcf
	else
	ln -s $file /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleName.amplicon.vcf
	fi


	bash /home/pipelines/master/shell/updatepipelineStatus.sh -q $queueID -s vep

	echo "running VEP"
	/home/pipelines/master/perl/bin/perl \
	/opt/vep/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \
	-i /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/$sampleName.amplicon.vcf \
	-o /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/$sampleName.amplicon.vep.vcf \
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

  bash /home/pipelines/master/shell/updatepipelineStatus.sh -q $queueID -s parseVEP

	echo "parse VEP results for new varView"
	python /home/pipelines/master/python/parseVEP.py \
	parseIllumina \
	-I /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/$sampleName.amplicon.vep.vcf \
	-o /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/$sampleName.amplicon.vep.parse.txt

	echo "filter VEP results"
	shopt -s nocasematch
	if [[ $sampleName =~ horizon ]]
	then
	awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 1 && $11 >=100) print}' /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/$sampleName.amplicon.vep.parse.txt \
	> /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/$sampleName.amplicon.vep.parse.filter.txt
	else
	awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 10 && $11 >=100) print}' /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/$sampleName.amplicon.vep.parse.txt \
	> /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/$sampleName.amplicon.vep.parse.filter.txt
	fi


	sed -i '1iGene\texon\tchr\tpos\tref\talt\tClassification\tType\tQuality\tAltVariantFreq\tRead Depth\tAltReadDepth\tConsequence\tSift\tPolyPhen\tHGVSc\tHGVSp\tdbSNPID\tpubmed' /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/$sampleName.amplicon.vep.parse.filter.txt

	##split amplicon file

	sampleCol=$(head -n 1 /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/amplicon.filtered.txt | awk -F"\t" -v sample=$sampleName '{for (i=1; i<= NF; i++) if($i==sample) print i}')

	if [ $sampleCol = "" ]
	then
	echo "amplicon info not found in amplicon file for sample:$sampleName"
	fi

	cut -f 1,$sampleCol /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/amplicon.filtered.txt > /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/$sampleName.amplicon.txt
	awk '$2 < 100' /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/$sampleName.amplicon.txt > /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/$sampleName.amplicon.lessThan100.txt


  bash /home/pipelines/master/shell/updatepipelineStatus.sh -q $queueID -s runCompleted

elif [ $instrumentID == "nextseq" ]
then

	runFolder=$(ls -d /home/$instrument/*_"$runID"_*)
	post=${runFolder##/home/$instrument/}
	runName=${post%%/Data/Intensities/BaseCalls/Alignment*}


	echo "runFolder is $runFolder"
	echo "post is $post"
	echo "runName is $runName"

	if [ ! -d /home/environments/$environmentID/"$instrumentID"Analysis/$runName ]
	then
		mkdir /home/environments/$environmentID/"$instrumentID"Analysis/$runName
	fi
	chmod 777 /home/environments/$environmentID/"$instrumentID"Analysis/$runName

	if [ ! -d /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID ]
	then
		mkdir /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID
	fi
	chmod 777 /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID

	##get runDate information##
	dateString=${runName%%_*}
	year=20${dateString:0:2}
	month=${dateString:2:2}
	day=${dateString:4:2}
	dateString=$year-$month-$day

	echo $dateString > /home/environments/$environmentID/"$instrumentID"Analysis/$runName/$sampleID/runDate.txt


	###Start processing file
	#for file in $runFolder/*.vcf
	for file in /home/environments/$environmentID/nextSeq_heme/$runName/$sampleID/*.vcf
	do
		echo $file
		if  [[ ${file} =~ $sampleID ]]
		then
			break
		fi
	done

	echo "testing $file"

	if [[ $file =~ .comb.vcf ]]
		##if the vcf is in the raw machine output form
	then
		echo "Processing $file "
		fileName=${file##*/}
		#sampleName=${fileName%%_*}
		sampleName=${fileName%%.comb*}

		echo "fileName is $fileName"
		echo "sampleName is $sampleName"


    bash /home/pipelines/master/shell/updatepipelineStatus.sh -q $queueID -s bedtools

			##filter vcf against ampliconRef
		if [ ! -z $ampliconRef ]
			then
			/opt/software/bedtools-2.17.0/bin/bedtools intersect -u -a $file -b $ampliconRef > /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleID/$sampleName.amplicon.vcf
		else
			ln -s $file /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleID/$sampleName.amplicon.vcf         #create symbolic link
		fi

    bash /home/pipelines/master/shell/updatepipelineStatus.sh -q $queueID -s vep

		echo "running VEP"
		/home/pipelines/master/perl/bin/perl \
		/opt/vep/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \
			-i /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleID/$sampleName.amplicon.vcf \
			-o /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleID/$sampleName.amplicon.vep.vcf \
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


			echo "parse VEP results for new varView"
			python /home/pipelines/master/python/parseVEP.py \
			parseIlluminaNextseq \
			-I /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleID/$sampleName.amplicon.vep.vcf \
			-o /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleID/$sampleName.amplicon.vep.parse.txt

			echo "filter VEP results"
			shopt -s nocasematch
			if [[ $sampleName =~ horizon ]]
				then
				awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 1 && $11 >=100) print}' /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleID/$sampleName.amplicon.vep.parse.txt \
				> /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleID/$sampleName.amplicon.vep.parse.filter.txt
			else
				awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 10 && $11 >=100) print}' /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleID/$sampleName.amplicon.vep.parse.txt \
				> /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleID/$sampleName.amplicon.vep.parse.filter.txt
			fi

			echo "running sed"
			sed -i '1iGene\texon\tchr\tpos\tref\talt\tClassification\tType\tQuality\tAltVariantFreq\tRead Depth\tAltReadDepth\tConsequence\tSift\tPolyPhen\tHGVSc\tHGVSp\tdbSNPID\tpubmed' /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleID/$sampleName.amplicon.vep.parse.filter.txt


      bash /home/pipelines/master/shell/updatepipelineStatus.sh -q $queueID -s samtools

			echo "generating "$sampleName" ampliconFile"
			/opt/samtools-1.4/samtools-1.4/samtools bedcov  /doc/ref/Heme/trusight-myeloid-amplicon-track.excludedNew.bed /home/environments/$environmentID/"$instrument"_heme/$runName/$sampleID/$sampleName.sort.bam | awk ' {print $4,"\t",int($13/($8-$7))} ' > /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleID/$sampleName.samtools.coverageDepth
			ampliconFile=$(ls /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.samtools.coverageDepth)
			###filter amplicon file
			if [ ! -z $excluded ]
			then
			grep -v -f $excluded $ampliconFile > /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleID/$sampleName.amplicon.filtered.txt
			else
			ln -s $ampliconFile /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleID/$sampleName.amplicon.filtered.txt
			fi



			##split amplicon file
			sampleCol=$(head -n 1 /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleID/$sampleName.amplicon.filtered.txt | awk -F"\t" '{print NF; exit}')
			echo "amplicon file has $sampleCol columns"
			if [ $sampleCol != "2" ]
				then
				echo "amplicon info not found in amplicon file for sample:$sampleName"
			fi


			cut -f 1,$sampleCol /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleID/$sampleName.amplicon.filtered.txt > /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleID/$sampleName.amplicon.txt
			awk '$2 < 100' /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleID/$sampleName.amplicon.txt > /home/environments/$environmentID/"$instrument"Analysis/$sampleID/$runName/$sampleName.amplicon.lessThan100.txt

			bash /home/pipelines/master/shell/updatepipelineStatus.sh -q $queueID -s runCompleted

	fi
	echo "Pipeline Finished!"


fi
echo "Pipeline Finished!"

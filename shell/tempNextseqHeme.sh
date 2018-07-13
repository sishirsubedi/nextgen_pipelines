if [ $# -eq 0 ]
	then
	echo "Usage: tempNextseqHeme.sh"
	echo "-r runID: four digits, last number of the run"
	echo "-a amplicon File in bed format (optional; required if -e is specified. Only variants within range specified in this file will be analyzed"
	echo "-e excluded amplicon file list, name only, one per line (optional)"
  echo "-n environment type ID"
	exit
fi

if test $# -gt 0
	then
	while getopts :i:r:a:e: opt
	do
	case $opt in
	r)
		runID=$OPTARG
		;;
	a)
		ampliconRef=$OPTARG
		;;
	e)
		excluded=$OPTARG
		;;
	n)
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

instrument=nextseq

if [ ! -d /home/environments/$environmentID/$instrument/*_"$runID"_*/ ]
	then
	echo "Error: result folder for run ID:$runID not found"
	exit
fi


if [ -z $instrument ]
	then
	echo "Error: instrument note defined"#must be specificed with -i option
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

if [ -z $runID ] || [ -z $instrument ]
	then
	echo "Error: -r runID and -i instrument are the required parameters"
	exit
fi


runFolder=$(ls -d /home/environments/$environmentID/$instrument/*_"$runID"_*)
post=${runFolder##/home/environments/$environmentID/$instrument/}
runName=${post%%/Data/Intensities/BaseCalls/Alignment*}


echo "runFolder is $runFolder"
echo "post is $post"
echo "runName is $runName"

if [ ! -d /home/environments/$environmentID/"$instrument"Analysis/$runName ]
	then
	mkdir /home/environments/$environmentID/"$instrument"Analysis/$runName
fi
chmod 777 /home/environments/$environmentID/"$instrument"Analysis/$runName

##get runDate information##
dateString=${runName%%_*}
year=20${dateString:0:2}
month=${dateString:2:2}
day=${dateString:4:2}
dateString=$year-$month-$day

echo $dateString > /home/environments/$environmentID/"$instrument"Analysis/$runName/runDate.txt





###Start processing file
#for file in $runFolder/*.vcf
for file in /home/environments/$environmentID/nextSeq_heme/$runName/*.vcf
	do
	echo "testing $file"
	#if [[ $file =~ .*_S[0-9]+\.vcf ]]
	if [[ $file =~ .comb.vcf ]]
		##if the vcf is in the raw machine output form
		then
		echo "Processing $file "
		fileName=${file##*/}
		#sampleName=${fileName%%_*}
		sampleName=${fileName%%.comb*}

		echo "fileName is $fileName"
		echo "sampleName is $sampleName"

		##filter vcf against ampliconRef
		if [ ! -z $ampliconRef ]
			then
			echo "running /opt/software/bedtools-2.17.0/bin/bedtools intersect -u -a $file -b $ampliconRef > /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.vcf"
			/opt/software/bedtools-2.17.0/bin/bedtools intersect -u -a $file -b $ampliconRef > /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.vcf
		else
			echo "running ln -s $file /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.vcf"
			ln -s $file /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.vcf         #create symbolic link
		fi

		echo "running VEP"
		/home/niyunyun/perl/bin/perl \
		/opt/vep/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \
		-i /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.vcf \
		-o /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.vep.vcf \
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
		-I /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.vep.vcf \
		-o /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.vep.parse.txt

		echo "filter VEP results"
		shopt -s nocasematch
		if [[ $sampleName =~ horizon ]]
			then
			awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 1 && $11 >=100) print}' /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.vep.parse.txt \
			> /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.vep.parse.filter.txt
		else
			awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 10 && $11 >=100) print}' /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.vep.parse.txt \
			> /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.vep.parse.filter.txt
		fi

		echo "running sed"
		sed -i '1iGene\texon\tchr\tpos\tref\talt\tClassification\tType\tQuality\tAltVariantFreq\tRead Depth\tAltReadDepth\tConsequence\tSift\tPolyPhen\tHGVSc\tHGVSp\tdbSNPID\tpubmed' /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.vep.parse.filter.txt



		echo "generating "$sampleName" ampliconFile"
		/opt/samtools-1.4/samtools-1.4/samtools bedcov  /doc/ref/Heme/trusight-myeloid-amplicon-track.excludedNew.bed /home/environments/$environmentID/nextSeq_heme/$runName/$sampleName.sort.bam | awk ' {print $4,"\t",int($13/($8-$7))} ' > /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.samtools.coverageDepth
		ampliconFile=$(ls /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.samtools.coverageDepth)
		###filter amplicon file
		if [ ! -z $excluded ]
		then
		grep -v -f $excluded $ampliconFile > /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.filtered.txt
		else
		ln -s $ampliconFile /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.filtered.txt
		fi



		##split amplicon file
		sampleCol=$(head -n 1 /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.filtered.txt | awk -F"\t" '{print NF; exit}')
		echo "amplicon file has $sampleCol columns"
		if [ $sampleCol != "2" ]
			then
			echo "amplicon info not found in amplicon file for sample:$sampleName"
		fi


		cut -f 1,$sampleCol /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.filtered.txt > /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.txt
		awk '$2 < 100' /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.txt > /home/environments/$environmentID/"$instrument"Analysis/$runName/$sampleName.amplicon.lessThan100.txt

	fi
done

echo "Pipeline Finished!"

if [ $# -eq 0 ]
	then
	echo "Usage: tempNextseqHeme.sh"
	echo "-r runID: four digits, last number of the run"
	echo "-a amplicon File in bed format (optional; required if -e is specified. Only variants within range specified in this file will be analyzed"
	echo "-e excluded amplicon file list, name only, one per line (optional)"
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

if [ ! -d /home/$instrument/*_"$runID"_*/ ]
	then
	echo "Error: result folder for run ID:$runID not found"
	exit
fi

#if [ ! -f /home/$instrument/*_"$runID"_*/Data/Intensities/BaseCalls/Alignment/AmpliconCoverage_M1.tsv ]
#	then
#	echo "Error: amplicon file for run ID:$runID not found"
#	exit
#fi

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

#ampliconFile=$(ls /home/$instrument/*_"$runID"_*/Data/Intensities/BaseCalls/Alignment/AmpliconCoverage_M1.tsv)
runFolder=$(ls -d /home/$instrument/*_"$runID"_*/)
post=${runFolder##/home/$instrument/}
runName=${post%%/Data/Intensities/BaseCalls/Alignment*}

echo "runFolder is $runFolder"
echo "post is $post"
echo "runName is $runName"

if [ ! -d /home/"$instrument"Analysis/$runName ]
	then
	mkdir /home/"$instrument"Analysis/$runName
fi
chmod 777 /home/"$instrument"Analysis/$runName

##get runDate information##
dateString=${runName%%_*}
year=20${dateString:0:2}
month=${dateString:2:2}
day=${dateString:4:2}
dateString=$year-$month-$day

echo $dateString > /home/"$instrument"Analysis/$runName/runDate.txt

###filter amplicon file
#if [ ! -z $excluded ]
#then
#grep -v -f $excluded $ampliconFile > /home/"$instrument"Analysis/$runName/amplicon.filtered.txt
#else
#ln -s $ampliconFile /home/"$instrument"Analysis/$runName/amplicon.filtered.txt
#fi

###Start processing file
#for file in $runFolder/*.vcf
for file in /home/nextSeq_heme/$runName/*.vcf
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
			echo "running /opt/software/bedtools-2.17.0/bin/bedtools intersect -u -a $file -b $ampliconRef > /home/"$instrument"Analysis/$runName/$sampleName.amplicon.vcf"
			/opt/software/bedtools-2.17.0/bin/bedtools intersect -u -a $file -b $ampliconRef > /home/"$instrument"Analysis/$runName/$sampleName.amplicon.vcf
		else
			echo "running ln -s $file /home/"$instrument"Analysis/$runName/$sampleName.amplicon.vcf"
			ln -s $file /home/"$instrument"Analysis/$runName/$sampleName.amplicon.vcf
		fi

		echo "running VEP"
		/home/niyunyun/perl/bin/perl \
		/opt/vep/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \
		-i /home/"$instrument"Analysis/$runName/$sampleName.amplicon.vcf \
		-o /home/"$instrument"Analysis/$runName/$sampleName.amplicon.vep.vcf \
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
		python /home/niyunyun/code/python/parseVEP.py \
		parseIlluminaNextseq \
		-I /home/"$instrument"Analysis/$runName/$sampleName.amplicon.vep.vcf \
		-o /home/"$instrument"Analysis/$runName/$sampleName.amplicon.vep.parse.txt

		echo "filter VEP results"
		shopt -s nocasematch
		if [[ $sampleName =~ horizon ]]
			then
			awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 1 && $11 >=100) print}' /home/"$instrument"Analysis/$runName/$sampleName.amplicon.vep.parse.txt \
			> /home/"$instrument"Analysis/$runName/$sampleName.amplicon.vep.parse.filter.txt
		else
			awk '{if(($7=="HIGH" || $7 =="MODERATE") && $10 >= 10 && $11 >=100) print}' /home/"$instrument"Analysis/$runName/$sampleName.amplicon.vep.parse.txt \
			> /home/"$instrument"Analysis/$runName/$sampleName.amplicon.vep.parse.filter.txt
		fi

		echo "running sed"
		sed -i '1iGene\texon\tchr\tpos\tref\talt\tClassification\tType\tQuality\tAltVariantFreq\tRead Depth\tAltReadDepth\tConsequence\tSift\tPolyPhen\tHGVSc\tHGVSp\tdbSNPID\tpubmed' /home/"$instrument"Analysis/$runName/$sampleName.amplicon.vep.parse.filter.txt
	
		##split amplicon file

		#sampleCol=$(head -n 1 /home/"$instrument"Analysis/$runName/amplicon.filtered.txt | awk -F"\t" -v sample=$sampleName '{for (i=1; i<= NF; i++) if($i==sample) print i}')
		#echo "sampleCol is $sampleCol"
		#if [ $sampleCol = "" ]
		#	then
		#	echo "amplicon info not found in amplicon file for sample:$sampleName"
		#fi
		
		#cut -f 1,$sampleCol /home/"$instrument"Analysis/$runName/amplicon.filtered.txt > /home/"$instrument"Analysis/$runName/$sampleName.amplicon.txt
		#awk '$2 < 100' /home/"$instrument"Analysis/$runName/$sampleName.amplicon.txt > /home/"$instrument"Analysis/$runName/$sampleName.amplicon.lessThan100.txt
		touch /home/"$instrument"Analysis/$runName/$sampleName.amplicon.txt
		touch /home/"$instrument"Analysis/$runName/$sampleName.amplicon.lessThan100.txt
		
	fi
done





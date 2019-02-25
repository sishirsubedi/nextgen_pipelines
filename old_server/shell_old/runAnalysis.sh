#!/bin/bash

##this is test line

if [ $# -eq 0 ]
then
	echo "-q queueID"
	echo "-e environment"
	echo "-u user"
	echo "-p password"
	exit
fi

if test $# -gt 0
	then
	while getopts :q:e:u:p: opt
	do
	case $opt in
  q)
		queueID=$OPTARG
		;;
	e)
	  environment=$OPTARG
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


 echo " $currentdate  INFO  - running runAnalysis for queueID - $queueID "

function updateStatus() {
user=$4
password=$5
database=$3
insertstatement="INSERT INTO pipelineStatus (queueID, plStatus, timeUpdated) VALUES ('$1','$2',now());"
mysql --user="$user" --password="$password" --database="$database" --execute="$insertstatement"
}


statement="select samples.sampleID, pipelineQueue.queueID, samples.runID, samples.sampleName, samples.coverageID, samples.callerID, assays.assayName, instruments.instrumentName, pipelineQueue.status from pipelineQueue join samples on samples.sampleID=pipelineQueue.sampleID join assays on assays.assayID = samples.assayID join instruments on instruments.instrumentID = samples.instrumentID where pipelineQueue.queueID='$queueID';"
variantstatement=""
ampliconstatement=""

## for checking successful database entry
variantfile_lines=-1
ampliconfile_lines=-1
current_sampleID=""

mysql --user="$user" --password="$password" --database="$environment" --execute="$statement" -N | while  read -r sampleID queueID  runID  sampleName coverageID callerID assay instrument status;
do

	current_sampleID=$sampleID

  if [ "$instrument" == "proton" ] || [ "$instrument" == "pgm" ]
  then
     echo " $currentdate   INFO  - getting Analysis -- instrument $instrument -- assay $assay -- run $runID -- coverage $coverageID -- caller $callerID  -- sample id is $sampleName"

     #enter not filtered variant data
     variantfile=$(ls /home/environments/$environment/"$instrument"Analysis/*$runID/$sampleName/"$callerID"/TSVC_variants.split.vep.parse.newVarView.txt)
		 #enter filtered variant data
		 #variantfile=$(ls /home/environments/$environment/"$instrumentID"Analysis/*$runID/$sampleID/"$callerID"/TSVC_variants.split.vep.parse.newVarView.filter.txt)

     variantfile_lines=$(wc -l < /home/environments/$environment/"$instrument"Analysis/*$runID/$sampleName/"$callerID"/TSVC_variants.split.vep.parse.newVarView.txt)

		 if [ ! -f $variantfile ]; then
		 		  updateStatus "$queueID" "ERROR:VariantFile" "$environment" "$user"  "$password"
          echo "Error Analsysis variant file not found !"
					exit

		 fi

		 variantstatement="load data local infile '$variantfile' into table sampleVariants (gene, exon, chr, pos, ref, alt, impact, type, quality, altFreq, readDepth, altReadDepth, consequence, Sift, PolyPhen, HGVSc, HGVSp, dbSNPID, pubmed) set sampleID = '$sampleID'"

		 runFolder=$(ls -d /home/environments/$environment/"$instrument"Analysis/*$runID)
		 runName=${runFolder##*/}

		 /opt/python3/bin/python3 /var/pipelines_"$environment"/python/loadSampleAmplicons.py \
		 													-i /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$coverageID"/amplicon.filter.txt \
															-o /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/"$coverageID"/amplicon.filter.v2.txt \
															-s $sampleID \
															-n $instrument\

    ampliconfile=$(ls /home/environments/$environment/"$instrument"Analysis/*$runID/$sampleName/"$coverageID"/amplicon.filter.v2.txt)

		ampliconfile_lines=$(wc -l < /home/environments/$environment/"$instrument"Analysis/*$runID/$sampleName/"$coverageID"/amplicon.filter.v2.txt)

		if [ ! -f $ampliconfile ]; then
		 		 updateStatus "$queueID" "ERROR:AmpliconFile" "$environment" "$user"  "$password"
		 		 echo "Error Analsysis amplicon file not found !"
		 		exit
		fi

		ampliconstatement="load data local infile '$ampliconfile' into table sampleAmplicons FIELDS TERMINATED BY '\t' (sampleID, gene, ampliconName,readDepth )"

  elif [ "$instrument" == "nextseq" ] || [ "$instrument" == "miseq" ]
  then


		 if [ "$instrument" == "nextseq" ]
		 then
     	echo " $currentdate    INFO  - getting Analysis  -- instrument $instrument -- assay $assay -- run $runID -- sample id is $sampleName"

		 	#enter not filtered variant data
     	variantfile=$(ls /home/environments/$environment/"$instrument"Analysis/*_"$runID"_*/$sampleName/variantAnalysis/$sampleName.amplicon.vep.parse.txt)

		 	#enter filtered variant data
     	#variantfile=$(ls /home/environments/$environment/"$instrumentID"Analysis/*_"$runID"_*/$sampleID/variantAnalysis/$sampleID.amplicon.vep.parse.filter.txt)

			variantfile_lines=$(wc -l < /home/environments/$environment/"$instrument"Analysis/*_"$runID"_*/$sampleName/variantAnalysis/$sampleName.amplicon.vep.parse.txt)

		 	if [ ! -f $variantfile ]; then
					 updateStatus "$queueID" "ERROR:VariantFile" "$environment" "$user"  "$password"
					 echo "Error Analsysis variant file not found !"
					exit
		 	fi

		 	variantstatement="load data local infile '$variantfile' into table sampleVariants (gene, exon, chr, pos, ref, alt, impact, type, quality, altFreq, readDepth, altReadDepth, consequence, Sift, PolyPhen, HGVSc, HGVSp, dbSNPID, pubmed) set sampleID = '$sampleID'"

     	#amplicon file
     	ampliconfile=$(ls /home/environments/$environment/"$instrument"Analysis/*_"$runID"_*/$sampleName/variantAnalysis/$sampleName.amplicon.lessThan100.txt)

			runFolder=$(ls -d /home/environments/$environment/"$instrument"Analysis/*_"$runID"_*)
 		  runName=${runFolder##*/}

			/opt/python3/bin/python3 /var/pipelines_"$environment"/python/loadSampleAmplicons.py \
															 -i /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.txt \
															 -o /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/variantAnalysis/$sampleName.amplicon.v2.txt \
															 -s $sampleID \
															 -n $instrument \

			#amplicon file
			ampliconfile=$(ls /home/environments/$environment/"$instrument"Analysis/*_"$runID"_*/$sampleName/variantAnalysis/$sampleName.amplicon.v2.txt)

			ampliconfile_lines=$(wc -l < /home/environments/$environment/"$instrument"Analysis/*_"$runID"_*/$sampleName/variantAnalysis/$sampleName.amplicon.v2.txt)

			if [ ! -f $ampliconfile ]; then
				updateStatus "$queueID" "ERROR:AmpliconFile" "$environment" "$user"  "$password"
				echo "Error Analsysis amplicon file not found !"
				exit
		  fi

      ampliconstatement="load data local infile '$ampliconfile' into table sampleAmplicons FIELDS TERMINATED BY '\t' (sampleID, gene, ampliconName,readDepth )"

	 elif  [ "$instrument" == "miseq" ]
	 then
		 echo " $currentdate  INFO  - getting Analysis  -- instrument $instrument -- assay $assay -- run $runID -- sample id is $sampleName"

      #enter variant without filtering data
      variantfile=$(ls /home/environments/$environment/"$instrument"Analysis/*_"$runID"_*/$sampleName/$sampleName.amplicon.vep.parse.txt)

			# #enter variant with filtering data
			# variantfile=$(ls /home/environments/$environment/"$instrument"Analysis/*_"$runID"_*/$sampleName/$sampleName.amplicon.vep.parse.filter.txt)

			variantfile_lines=$(wc -l < /home/environments/$environment/"$instrument"Analysis/*_"$runID"_*/$sampleName/$sampleName.amplicon.vep.parse.txt)

 		 if [ ! -f $variantfile ]; then
 					 updateStatus "$queueID" "ERROR:VariantFile" "$environment" "$user"  "$password"
 					 echo "Error Analsysis variant file not found !"
 					exit
 		 fi

 		 variantstatement="load data local infile '$variantfile' into table sampleVariants  (gene, exon, chr, pos, ref, alt, impact, type, quality, altFreq, readDepth, altReadDepth, consequence, Sift, PolyPhen, HGVSc, HGVSp, dbSNPID, pubmed) set sampleID = '$sampleID'"

		 runFolder=$(ls -d /home/environments/$environment/"$instrument"Analysis/*_"$runID"_*)
		 runName=${runFolder##*/}

			/opt/python3/bin/python3 /var/pipelines_"$environment"/python/loadSampleAmplicons.py \
															 -i /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.txt \
															 -o /home/environments/$environment/"$instrument"Analysis/$runName/$sampleName/$sampleName.amplicon.v2.txt \
															 -s $sampleID \
															 -n $instrument \

	     #amplicon file
			 ampliconfile=$(ls /home/environments/$environment/"$instrument"Analysis/*_"$runID"_*/$sampleName/$sampleName.amplicon.v2.txt)

			 ampliconfile_lines=$(wc -l < /home/environments/$environment/"$instrument"Analysis/*_"$runID"_*/$sampleName/$sampleName.amplicon.v2.txt)

		   if [ ! -f $ampliconfile ]; then
				  updateStatus "$queueID" "ERROR:AmpliconFile" "$environment" "$user"  "$password"
			    echo "Error Analsysis amplicon file not found !"
					exit
			fi

		  ampliconstatement="load data local infile '$ampliconfile' into table sampleAmplicons FIELDS TERMINATED BY '\t' (sampleID, gene, ampliconName,readDepth )"

	 fi


  fi

mysql --user="$user" --password="$password" --database="$environment" --execute="$variantstatement"
mysql --user="$user" --password="$password" --database="$environment" --execute="$ampliconstatement"


### check data entry into database

variants_check_statement="select count(*) from sampleVariants where sampleID=$current_sampleID"
amplicons_check_statement="select count(*) from sampleAmplicons where sampleID=$current_sampleID"

variantdb_lines=$(mysql --user="$user" --password="$password" --database="$environment" --execute="$variants_check_statement" -N)
amplicondb_lines=$(mysql --user="$user" --password="$password" --database="$environment" --execute="$amplicons_check_statement" -N)

echo "variant file lines is "$variantfile_lines" and amplicon file line is "$ampliconfile_lines" "
echo "variant db lines is "$variantdb_lines" and amplicon db line is "$amplicondb_lines" "

if [ "$variantfile_lines" == "$variantdb_lines" ]; then
		echo "sample variants entered in database successfully."
else
	echo "ERROR: sample variants entered in database FAILED."
	updateStatus "$queueID" "ERROR:VariantFile" "$environment" "$user"  "$password"
fi

if [ "$ampliconfile_lines" == "$amplicondb_lines" ]; then
		echo "sample amplicons entered in database successfully."
else
	echo "ERROR: sample amplicons entered in database FAILED."
	updateStatus "$queueID" "ERROR:AmpliconFile" "$environment" "$user"  "$password"
fi

#### get pipeline entering user
usergetstatement="select email from users where userID in (select enteredBy from samples as t1  join pipelineQueue as t2 on  t1.sampleID=t2.sampleID where t2.queueID ='$queueID');"
useremail=$(mysql --user="$user" --password="$password" --database="$environment" --execute="$usergetstatement" -N)

echo " $currentdate   INFO  - Emailing $useremail with completed message"

message=$'Pipeline Completed for-
				assay : '$assay' ,
				instrument : '$instrument' ,
				runID : '$runID' ,
				sampleName : '$sampleName' ,
				coverageID : '$coverageID' ,
				callerID : '$callerID' ,
				environment : '$environment' ,
				queueID : '$queueID' '
echo $message

echo $message |  mailx  -s "Pipeline completed for sample:$sampleName , run:$runID" "$useremail"

done

updateStatus "$queueID" "pipelineCompleted" "$environment" "$user"  "$password"

/opt/python3/bin/python3 /var/pipelines_"$environment"/python/logPipelineRun.py -q "$queueID" -u "$user" -p "$password" -d "$environment"

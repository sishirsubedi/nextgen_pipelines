#!/bin/bash

##this is test line

if [ $# -eq 0 ]
then
	echo "-q queueID"
	echo "-e invironmentID"
	exit
fi

if test $# -gt 0
	then
	while getopts :q:e: opt
	do
	case $opt in
  q)
		queueID=$OPTARG
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


###database credentials --
user=hhadmin
password=ngs3127
database="$environmentID"


function updateStatus() {

user=hhadmin
password=ngs3127
database=test
insertstatement="INSERT INTO pipelineStatus (queueID, plStatus, timeUpdated) VALUES ('$1','$2',now());"
mysql --user="$user" --password="$password" --database="$database" --execute="$insertstatement"
}



statement="select ID, queueID,  runID,  sampleID, coverageID, vcallerID,  assayID, instrumentID, environmentID, status from sampleAnalysisQueue where queueID='$queueID' order by queueID asc limit 1;"

variantstatement=""
coveragestatement=""
ampliconstatement=""

mysql --user="$user" --password="$password" --database="$database" --execute="$statement" -N | while  read -r ID queueID  runID  sampleID coverageID vcallerID assayID instrumentID environmentID status;
do
  echo " $ID : $queueID  :$runID  :$sampleID :$coverageID :$vcallerID :$assayID :$instrumentID :$environmentID: $status"

  if [ "$instrumentID" == "proton" ] || [ "$instrumentID" == "pgm" ]
  then

     echo "getting Analysis -- instrument $instrumentID -- assay $assayID -- run $runID -- sample id is $sampleID"

     #enter variant data
     variantfile=$(ls /home/environments/$environmentID/"$instrumentID"Analysis/*$runID/$sampleID/variantCaller_out."$vcallerID"/TSVC_variants.split.vep.parse.newVarView.filter.txt)

		 if [ ! -f $variantfile ]; then
		 		  updateStatus "$queueID" "ERROR:VariantFile" "$environmentID"
          echo "Error Analsysis variant file not found !"
					exit
		 fi

		 variantstatement="load data local infile '$variantfile' into table data ignore 1 lines (gene, exons, chr, pos, ref, alt, genotype, type, quality, altFreq, readDP, altReadDP, Consequence, sift, PolyPhen, HGVSc, HGVSp, dbSNPID, pubmed) set sampleID = '$ID', assay = '$assayID'"

     #amplicon filebash
     ampliconfile=$(ls /home/environments/$environmentID/"$instrumentID"Analysis/*$runID/$sampleID/coverageAnalysis_out."$coverageID"/amplicon.lessThan100.txt)

		 if [ ! -f $ampliconfile ]; then
					 updateStatus "$queueID" "ERROR:AmpliconFile" "$environmentID"
					 echo "Error Analsysis amplicon file not found !"
					exit
		 fi

		 coveragestatement="load data local infile '$ampliconfile' into table amplicon (ampliconName, ampliconCov) set sampleID = '$ID', assay = '$assayID'"

     #total amplicon number
     totalAmpliconCount=$(wc -l  /home/environments/$environmentID/"$instrumentID"Analysis/*$runID/$sampleID/coverageAnalysis_out."$coverageID"/amplicon.filter.txt | cut -d ' ' -f 1 )
     #total failed amplicon number
     failedAmpliconCount=$(wc -l /home/environments/$environmentID/"$instrumentID"Analysis/*$runID/$sampleID/coverageAnalysis_out."$coverageID"/amplicon.lessThan100.txt | cut -d ' ' -f 1)
     ampliconstatement="insert into ampliconCount values ('$ID', '$totalAmpliconCount', '$failedAmpliconCount');"



  elif [ "$instrumentID" == "nextseq" ] || [ "$instrumentID" == "miseq" ]
  then

     echo "getting Analysis  -- instrument $instrumentID -- assay $assayID -- run $runID -- sample id is $sampleID"

     #enter variant data
     variantfile=$(ls /home/environments/$environmentID/"$instrumentID"Analysis/*_"$runID"_*/$sampleID/$sampleID.amplicon.vep.parse.filter.txt)

		 if [ ! -f $variantfile ]; then
					 updateStatus "$queueID" "ERROR:VariantFile" "$environmentID"
					 echo "Error Analsysis variant file not found !"
					exit
		 fi

		 variantstatement="load data local infile '$variantfile' into table data ignore 1 lines (gene, exons, chr, pos, ref, alt, genotype, type, quality, altFreq, readDP, altReadDP, Consequence, sift, PolyPhen, HGVSc, HGVSp, dbSNPID, pubmed) set sampleID = '$ID', assay = '$assayID'"

     #amplicon file
     ampliconfile=$(ls /home/environments/$environmentID/"$instrumentID"Analysis/*_"$runID"_*/$sampleID/$sampleID.amplicon.lessThan100.txt)

		 if [ ! -f $ampliconfile ]; then
					 updateStatus "$queueID" "ERROR:AmpliconFile" "$environmentID"
					 echo "Error Analsysis amplicon file not found !"
					exit
		 fi

		 coveragestatement="load data local infile '$ampliconfile' into table amplicon (ampliconName, ampliconCov) set sampleID = '$ID', assay = '$assayID'"

     #total amplicon number
     totalAmpliconCount=$(wc -l  /home/environments/$environmentID/"$instrumentID"Analysis/*_"$runID"_*/$sampleID/$sampleID.amplicon.txt | cut -d ' ' -f 1 )
     #total failed amplicon number
     failedAmpliconCount=$(wc -l /home/environments/$environmentID/"$instrumentID"Analysis/*_"$runID"_*/$sampleID/$sampleID.amplicon.lessThan100.txt | cut -d ' ' -f 1)
     ampliconstatement="insert into ampliconCount values ('$ID', '$totalAmpliconCount', '$failedAmpliconCount');"

  fi

mysql --user="$user" --password="$password" --database="$database" --execute="$variantstatement"
mysql --user="$user" --password="$password" --database="$database" --execute="$coveragestatement"
mysql --user="$user" --password="$password" --database="$database" --execute="$ampliconstatement"


#### get pipeline entering user

usergetstatement="select email from networkIDTOemail where networkID in (select enteredBy from Samples as t1  join sampleAnalysisQueue as t2 on  t1.ID=t2.ID where t2.queueID ='$queueID');"
useremail=$(mysql --user="$user" --password="$password" --database="$database" --execute="$usergetstatement" -N)

echo "Emailing $useremail with completed message"

message=$'Pipeline Completed for-
				assayID : '$assayID' ,
				instrumentID : '$instrumentID' ,
				runID : '$runID' ,
				sampleID : '$sampleID' ,
				coverageID : '$coverageID' ,
				callerID : '$callerID' ,
				environmentID : '$environmentID' ,
				queueID : '$queueID' '
echo $message
mail -s "Pipeline completed for sample:$sampleID , run:$runID"   "$useremail" << MSG_BODY_HERE
Pipeline Completed for-
				assayID : $assayID ,
				instrumentID : $instrumentID ,
				runID : $runID ,
				sampleID : $sampleID ,
				coverageID : $coverageID ,
				callerID : $vcallerID ,
				environmentID : $environmentID ,
				queueID : $queueID
MSG_BODY_HERE

done

updateStatus "$queueID" "pipelineCompleted" "$environmentID"

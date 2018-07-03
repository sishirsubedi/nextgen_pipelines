#!/bin/bash

if [ $# -eq 0 ]
then
	echo "-q queueID"
	exit
fi

if test $# -gt 0
	then
	while getopts :q: opt
	do
	case $opt in
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

user=root
password=molSeq3127
database=test
table=sampleAnalysisQueue

statement="select ID, queueID,  runID,  sampleID, coverageID, vcallerID,  assayID, instrumentID, environmentID, status from $table where queueID='$queueID' order by queueID asc limit 1;"

variantstatement=""
coveragestatement=""
ampliconstatement=""

mysql --user="$user" --password="$password" --database="$database" --execute="$statement" -N | while  read -r ID queueID  runID  sampleID coverageID vcallerID assayID instrumentID environmentID status;
do
  echo " $ID : $queueID  :$runID  :$sampleID :$coverageID :$vcallerID :$assayID :$instrumentID :$environmentID: $status"

  if [ "$instrumentID" == "proton" ] || [ "$instrumentID" == "pgm" ]
  then

     echo "proton or pgm -- sample id is $sampleID"

     #enter variant data
     variantfile=$(ls /home/environments/$environmentID/"$instrumentID"Analysis/*$runID/variantCaller_out."$vcallerID"/$sampleID/TSVC_variants.split.vep.parse.newVarView.filter.txt)
     variantstatement="load data local infile '$variantfile' into table data ignore 1 lines (gene, exons, chr, pos, ref, alt, genotype, type, quality, altFreq, readDP, altReadDP, Consequence, sift, PolyPhen, HGVSc, HGVSp, dbSNPID, pubmed) set sampleID = '$ID', assay = '$assayID'"

     #amplicon filebash
     ampliconfile=$(ls /home/environments/$environmentID/"$instrumentID"Analysis/*$runID/coverageAnalysis_out."$coverageID"/$sampleID/amplicon.lessThan100.txt)
     coveragestatement="load data local infile '$ampliconfile' into table amplicon (ampliconName, ampliconCov) set sampleID = '$ID', assay = '$assayID'"

     #total amplicon number
     totalAmpliconCount=$(wc -l  /home/environments/$environmentID/"$instrumentID"Analysis/*$runID/coverageAnalysis_out."$coverageID"/$sampleID/amplicon.filter.txt | cut -d ' ' -f 1 )
     #total failed amplicon number
     failedAmpliconCount=$(wc -l /home/environments/$environmentID/"$instrumentID"Analysis/*$runID/coverageAnalysis_out."$coverageID"/$sampleID/amplicon.lessThan100.txt | cut -d ' ' -f 1)
     ampliconstatement="insert into ampliconCount values ('$ID', '$totalAmpliconCount', '$failedAmpliconCount');"



  elif [ "$instrumentID" == "nextseq" ] || [ "$instrumentID" == "miseq" ]
  then

     echo "nextseq or miseq"

     #enter variant data
     variantfile=$(ls /home/environments/$environmentID/"$instrumentID"Analysis/*_"$runID"_*/$sampleID.amplicon.vep.parse.filter.txt)
     variantstatement="load data local infile '$variantfile' into table data ignore 1 lines (gene, exons, chr, pos, ref, alt, genotype, type, quality, altFreq, readDP, altReadDP, Consequence, sift, PolyPhen, HGVSc, HGVSp, dbSNPID, pubmed) set sampleID = '$ID', assay = '$assayID'"

     #amplicon file
     ampliconfile=$(ls /home/environments/$environmentID/"$instrumentID"Analysis/*_"$runID"_*/$sampleID.amplicon.lessThan100.txt)
     coveragestatement="load data local infile '$ampliconfile' into table amplicon (ampliconName, ampliconCov) set sampleID = '$ID', assay = '$assayID'"

     #total amplicon number
     totalAmpliconCount=$(wc -l  /home/environments/$environmentID/"$instrumentID"Analysis/*_"$runID"_*/$sampleID.amplicon.txt | cut -d ' ' -f 1 )
     #total failed amplicon number
     failedAmpliconCount=$(wc -l /home/environments/$environmentID/"$instrumentID"Analysis/*_"$runID"_*/$sampleID.amplicon.lessThan100.txt | cut -d ' ' -f 1)
     ampliconstatement="insert into ampliconCount values ('$ID', '$totalAmpliconCount', '$failedAmpliconCount');"

  fi

mysql --user="$user" --password="$password" --database="$database" --execute="$variantstatement"
mysql --user="$user" --password="$password" --database="$database" --execute="$coveragestatement"
mysql --user="$user" --password="$password" --database="$database" --execute="$ampliconstatement"

done

bash /home/pipelines/master/shell/updatepipelineStatus.sh -q $queueID -s pipelineCompleted

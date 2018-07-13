runID=$1
sampleID=$2

query="select uid from molseq.hemeSample where instrument='miseq' and runID='$runID' and sampleID='$sampleID'"

uid=$(mysql molseq -u webuser --password=molSeq3127 -se "$query")

if [[ $uid == "" ]]
then
echo "Aborted: sample info does not exist."
exit
fi

###check if sample ID corresponds to sampleID in amplicon info file###
if [ -f /home/miseq/*_"$runID"_*/Data/Intensities/BaseCalls/Alignment/AmpliconCoverage_M1.filtered.tsv ]
then
ampliconFile=$(ls /home/miseq/*_"$runID"_*/Data/Intensities/BaseCalls/Alignment/AmpliconCoverage_M1.filtered.tsv)
elif [ -f /home/miseq/*_"$runID"_*/Data/Intensities/BaseCalls/Alignment/AmpliconCoverage_M1.tsv ]
then
ampliconFile=$(ls /home/miseq/*_"$runID"_*/Data/Intensities/BaseCalls/Alignment/AmpliconCoverage_M1.tsv)
else
echo "amplicon file not found for runID:$runID"
echo "amplicon info not entered for runID:$runID, sampleID:$sampleID"
exit
fi

sampleIDSimple=${sampleID%_*}
sampleCol=$(head -n 1 $ampliconFile | awk -F"\t" -v sample=$sampleIDSimple '{for (i=1; i<= NF; i++) if($i==sample) print i}')

if [ $sampleCol = "" ]
then
echo "amplicon info not found in amplicon file for sample:$sampleID"
echo "amplicon info not entered for runID:$runID, sampleID:$sampleID"
exit
fi

#sampleCol=$((sampleCol+1))

cut -f 1,$sampleCol $ampliconFile > /home/scratch/amplicon.temp.txt

cmd="load data local infile '/home/scratch/amplicon.temp.txt' into table molseq.hemeAmplicon ignore 1 lines (ampliconName, ampliconCov) set sampleID=$uid"

mysql -e "$cmd" -u webuser --password=molSeq3127
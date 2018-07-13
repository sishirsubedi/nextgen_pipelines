runID=$1
ampliconFile=$2
outFile=$3

colNum=$(awk '{print NF}' $ampliconFile |head -n 1) 

sampleNum=$((colNum-1))

ampliconDir=${ampliconFile%/*}



n=2
if [ $colNum != "" ]
then
while [ $n -le $colNum ]
do
sampleName=$(cut -f $n $ampliconFile |head -n 1)
statFile=$(ls $ampliconDir/$sampleName*.stats)
mappedReads=$(sed -n '3p' $statFile |cut -d ' ' -f 1)
failed=$(tail -n +2 $ampliconFile |awk -F'\t' -v OFS='\t' -v n=$n '{if($n<100) print}' |wc -l)
echo -e "$runID\t$sampleName\t$sampleNum\t$mappedReads\t$failed" >> $outFile
n=$((n+1))
done
fi
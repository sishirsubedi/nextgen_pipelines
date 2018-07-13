dir=$1
currentDir=$(pwd)

## not enough arguments##
if [ $# -lt 3 ]
then
echo "Not enough arguments"
echo "Usage: sh /home/MDXLAB/blastMiSeq.sh runID BlastType SampleID......"
echo "runID: first field in Run folder name"
echo "BlastType: 50K, 500K, control, full"
echo "SampleID: The sample name typed into your SampleSheet" 
exit
fi

##not the correct folder name##
if [ ! -d /home/miseq/"$dir"_*/ ]
then
echo "Error: run folder starting with $dir not found; nothing was done"
exit
fi

blastType=$(echo $2 |awk '{print tolower($0)}')
#determine blast type##
if [[ $blastType == "50k" ]]
then
workDir="/home/MDXLAB/Blast50K"
script="Blast50K.pl"
echo "performing blast 50K"
elif [[ $blastType == "500k" ]]
then
workDir="/home/MDXLAB/Blast500K"
script="Blast500K.pl"
echo "performing blast 500K"
elif [[ $blastType == "control" ]] || [[ $blastType == "controls" ]]
then
script="BlastControls.pl"
workDir="/home/MDXLAB/BlastControls"
echo "performing blast control"
elif [[ $blastType == "full" ]]
then
script="BlastFull.pl"
workDir="/home/MDXLAB/BlastFull"
echo "performing blast full"
else
echo "Error: $blastType is not a recognized blast type, please choose between: 50k, 500k, control, full"
exit
fi

##link fastq files to the corresponding blast dir##
fromDir=$(ls -d /home/miseq/"$dir"_*/Data/Intensities/BaseCalls)
n=3
while [[ $n -le $# ]]
do
sample=$(eval "echo \$$n")
if [ ! -f $fromDir/$sample* ]
then
echo "results for $sample do not exist, please double check"
echo $sample 
exit
else
ln -s $fromDir/$sample*.fastq.gz $workDir/
fi
(( n++ ))
done > $blastType.txt

##extract and perform blast##
cd $workDir

for file in *.fastq.gz
do
fileName=${file%.gz}
#gunzip -c $file > $fileName
done

#perl $script

##move files back to current dir##
cd $currentDir
while read p
do
mv $workDir/"$p"_* ./
done < $blastType.txt


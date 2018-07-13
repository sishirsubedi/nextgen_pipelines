dir=$1

if [ ! -d /home/miseq/"$dir"_*/ ]
then
echo "Error: run folder starting with $dir not found; nothing copied"
exit
fi

fromDir=$(ls -d /home/miseq/"$dir"_*)
ls $fromDir/Data/Intensities/BaseCalls/*.fastq.gz > fileList.txt

grep -v "Blank" fileList.txt |grep -v "Undetermined" > fileListSample.txt

echo "copying fastq.gz files from $fromDir to current directory and extracting them"

while read p
do
file=${p##*/}
cp $p ./
gunzip $file
done < fileListSample.txt

rm fileList.txt
rm fileListSample.txt
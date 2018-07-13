#this script runs bcl2fastq conversion and concatenate the R1 and R2 fastq files
#run bcl2fastq in the current directory

bcl2fastq --runfolder-dir ./ --output-dir ./out1 --no-lane-splitting

#unzip resulting R1 fastq files
mkdir ./out2
cd out1
#get a list of sample names
ls *.gz|cut -d '_' -f 1 |uniq > ../sampleList.txt
cd ..
#unzipping
while read p
do
echo unzipping $p
gunzip -c out1/$p*R1*.fastq.gz > out2/$p.fastq
done < sampleList.txt 
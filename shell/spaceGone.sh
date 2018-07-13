#! /bin/bash
#####take a directory and replace space with _ in all sub-directory and files in it#######
current_time=$(date +"%m-%d-%Y-%H-%M")
dir=$1
#rename directories first
find $dir -name "* *" -type d  > $current_time.D.txt
while read p
do
newDir=${p// /_}
mv "$p" $newDir
done < $current_time.D.txt
#rename files second
find $dir -name "* *" -type f > $current_time.F.txt
while read p
do
newFile=${p// /_}
mv "$p" $newFile
done < $current_time.F.txt
#delete temp files
#rm $current_time.D.txt
#rm $current_time.F.txt
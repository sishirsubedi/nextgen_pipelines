file=$1
file_name=${file##*/}
dir=${file%/*}
echo $file_name
echo $dir
if [ $file_name == "RunCompletionStatus.xml" ]
then
if [ -f $dir/SampleSheet.csv ]
then
#run bcl2fastq
ulimit -n 500000
echo "Hard limit" >> $dir/fileLimit.txt
ulimit -Hn >> $dir/fileLimit.txt
echo "Soft limit" >> $dir/fileLimit.txt
ulimit -Sn >> $dir/fileLimit.txt
/usr/local/bin/bcl2fastq --runfolder-dir $dir --output-dir $dir/out1 --no-lane-splitting > $dir/bcl2fastq.log 2>&1
#count read number

for file in $dir/out1/*.fastq.gz
do
sample=${file##*/}
echo -ne $sample"\t" >> $dir/out1/readCounts.txt
gunzip -c $file |wc -l >> $dir/out1/readCounts.txt
done

awk '{print $1, $2/4}' $dir/out1/readCounts.txt > $dir/out1/readCounts.final.txt
fi

#after the run, remove the dir from the monitor list
sed -i "\:$dir IN_CREATE sh /home/niyunyun/code/shell/bcl2fastqCron.sh \$@/\$#:d" /var/spool/incron/root
fi


#exec >> $dir/postProcess.log 2>&1
file=$1
file_name=${file##*/}
dir=${file%/*}
#echo $file_name
#echo $dir
if [ $file_name == "CompletedJobInfo.xml" ]
then
#run filter vcf
sleep 15m
echo "Start post processing ..." >> $dir/Data/Intensities/BaseCalls/Alignment/postProcess.log
sampleNum=$(ls $dir/Data/Intensities/BaseCalls/Alignment/*.vcf |wc -l)
echo "Total samples to process: "$sampleNum >> $dir/Data/Intensities/BaseCalls/Alignment/postProcess.log
if [[ $sampleNum -ge 1 ]]
then
for file in $dir/Data/Intensities/BaseCalls/Alignment/*.vcf
do
fileName=${file##*/}
sample=${fileName%.vcf}
grep "#CHROM" $file > $dir/Data/Intensities/BaseCalls/Alignment/$sample.header
/opt/software/bedtools-2.17.0/bin/bedtools intersect -u -a $file -b /home/doc/ref/Heme/trusight-myeloid-amplicon-track.excluded.bed > $dir/Data/Intensities/BaseCalls/Alignment/$sample.amp.vcf
cat $dir/Data/Intensities/BaseCalls/Alignment/$sample.header $dir/Data/Intensities/BaseCalls/Alignment/$sample.amp.vcf > $dir/Data/Intensities/BaseCalls/Alignment/$sample.amplicon.vcf
rm $dir/Data/Intensities/BaseCalls/Alignment/$sample.amp.vcf
rm $dir/Data/Intensities/BaseCalls/Alignment/$sample.header
done
fi

grep -v -f /home/doc/ref/Heme/excludedAmplicons.txt $dir/Data/Intensities/BaseCalls/Alignment/AmpliconCoverage_M1.tsv > $dir/Data/Intensities/BaseCalls/Alignment/AmpliconCoverage_M1.filtered.tsv

#after the run, remove the dir from the monitor list
sed -i "\:$dir IN_CREATE sh /home/niyunyun/code/shell/filterHemeVcf.sh \$@/\$#:d" /var/spool/incron/root
fi


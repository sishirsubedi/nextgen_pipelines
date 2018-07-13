#!/bin/bash
if [ $# -eq 0 ]
then
echo "Usage: hemeCosmic.sh runID"
echo "runID is the third field of the run folder name"
exit
fi

runID=$1
dir=$(ls -d /home/vmUser/miseq_out/*_*\_$runID\_*-*)
for file in $dir/*.tsv
do
fileName=${file##*/}
sample=${fileName%.tsv}
echo "processing "$sample
head -n 1 $file |awk '{ sub(/\r$/,""); print }' | awk -F'\t' -v OFS='\t' '{print $1,"chr","pos","ref","alt",$5,$6,$7,$10,$13,$14,$15,$22,$31,$32,$34,$35,$36,$37,$39,$40,$41, $42,$43,$44,$49,$50,$51,$55,$56,$57,$58,$59, "cosmicV72"}'> $dir/$sample.header
join -t$'\t' -a 1 -1 1 -2 1 \
<(awk '{ sub(/\r$/,""); print }' $file | awk -F'\t' -v OFS='\t' '{split($2, a, ">"); split($2, b, "/"); print "chr"$3"|"$4"|"a[1]"|"b[2], $1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45,$46,$47,$48,$49,$50,$51,$52,$53,$54,$55,$56,$57,$58,$59,$60,$61,$62,$63,$64,$65,$66,$67,$68,$69,$70,$71}' |tail -n +2 |sort -k 1,1) \
<(cut -f 1,2 /home/doc/ref/Heme/CosmicCodingMutsV72.inHeme.sort.txt) \
|awk -F'\t' -v OFS='\t' '{split($1, a, "|"); print ($2, a[1], a[2], a[3], a[4], $3,$4,$5,$8,$11,$12,$13,$20,$29,$30,$32,$33,$34,$35,$37, $38,$39,$40,$41,$42,$47,$48,$49,$53,$54,$55,$56,$57,$70)}' >  $dir/$sample.newCosmic.tsv
cat $dir/$sample.header $dir/$sample.newCosmic.tsv > $dir/$sample.newCosmic.final.tsv
rm $dir/$sample.newCosmic.tsv
rm $dir/$sample.header
done

miseqDir=$(ls -d /home/miseq/*_*\_$runID\_*-*/Data/Intensities/BaseCalls/Alignment/)
col=$(awk '{print NF}' $miseqDir/AmpliconCoverage_M1.tsv |head -n 1)

n=2

while [[ $n -le $col ]]
do
sample=$(head -n 1 $miseqDir/AmpliconCoverage_M1.tsv |cut -f $n)
awk -v OFS='\t' -v n=$n '{if($n < 100) print $1}' $miseqDir/AmpliconCoverage_M1.tsv > $dir/$sample.amplicon.txt
sed -i '1iInsufficiently Covered Amplicons' $dir/$sample.amplicon.txt
sed -i "\:NRAS.exon.3.line.6.chr1.115256420.115256599_tile_1:d" $dir/$sample.amplicon.txt
sed -i "\:DNMT3A.CDS.6.line.68.chr2.25475062.25475066_tile_1:d" $dir/$sample.amplicon.txt
sed -i "\:CUX1.CDS.32.line.153.chr7.101459311.101459373_tile_1:d" $dir/$sample.amplicon.txt
n=$((n+1))
done

for file in $dir/*.newCosmic.final.tsv
do
fileName=${file##*/}
sample=${fileName%%.*}
cat $file $dir/$sample*amplicon.txt > $dir/$sample.final.txt
rm $dir/$sample.amplicon.txt
mv $file $dir/$sample.newCosmic.txt
done

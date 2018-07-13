#save current time for log use
current_time=$(date +"%m-%d-%Y-%H-%M")

###find all .bcl, .bcl.bgzf, .bcl.bgzf.bci, .locs that are older than 90 days in /home/nextseq and /home/miseq folders###
find /home/miseq -mtime +90 -name *.bcl >> /home/niyunyun/cronLog/$current_time.file.txt
find /home/miseq -mtime +90 -name *.locs >> /home/niyunyun/cronLog/$current_time.file.txt
find /home/miseq -mtime +90 -name *.tif >> /home/niyunyun/cronLog/$current_time.file.txt
find /home/nextseq -mtime +90 -name *.locs >> /home/niyunyun/cronLog/$current_time.file.txt
find /home/nextseq -mtime +90 -name *.bcl.bgzf >> /home/niyunyun/cronLog/$current_time.file.txt
find /home/nextseq -mtime +90 -name *.bcl.bgzf.bci >> /home/niyunyun/cronLog/$current_time.file.txt


###delete the above files###
while read p
do
rm $p
done < /home/niyunyun/cronLog/$current_time.file.txt > /home/niyunyun/cronLog/$current_time.delete.log 2>&1

####change file permission of /home/nextseq and /home/miseq ####
chmod -R 775 /home/nextseq
chmod -R 775 /home/miseq

####rsync the new miseq directory structure to vmUser's home dir#####
#rsync -dv /home/miseq/ /home/vmUser/miseq_out > /home/niyunyun/cronLog/$current_time.rsync.log 2>&1
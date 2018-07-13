dir=$1
if [ -d $dir ]
then
echo $dir IN_CREATE sh /home/niyunyun/code/shell/bcl2fastqCron.sh '$@/$#' >> /var/spool/incron/root
fi

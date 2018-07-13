dir=$1
if [ -d $dir ]
then
#echo $dir IN_CREATE sh /home/niyunyun/code/shell/filterHemeVcf.sh '$@/$#' >> /var/spool/incron/root
echo $dir IN_CREATE bash /home/niyunyun/code/shell/hemePipeline.sh '$@/$#' >> /var/spool/incron/root
fi

#local back up
currentTime=$(date +%b-%d-%Y-%H-%M)
rsync -avrP /var/www/html /home/localBak/ > /home/niyunyun/cronlog/$currentTime.localBak.log 2>&1
rsync -avrP /home/niyunyun/code /home/localBak/ >> /home/niyunyun/cronlog/$currentTime.localBak.log 2>&1

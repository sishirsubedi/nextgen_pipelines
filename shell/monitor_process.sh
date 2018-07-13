#!/bin/bash
#!/bin/bash

if [ $# -eq 0 ]
then
	echo "Usage: monitor_process.sh"
	echo "-p processID: provide process ID to monitor usuage"
	exit
fi

processID=$1

echo "monitoring CPU_MEMORY usuage for process ID $processID "
filename="${processID}_Monitoring.txt"
touch $filename
i=0
# while [ $i -le 10 ] #runs 10 minutes
while [ ! ps -p $processID > /dev/null 2>&1; echo "..."]
do
    i=$((i+1))
    echo "monitoring process $processID ... updated $i times"
    sleep 10s
    ps -p $processID -o %cpu,%mem,cmd >> $filename
    echo -e '\n'>> $filename
done

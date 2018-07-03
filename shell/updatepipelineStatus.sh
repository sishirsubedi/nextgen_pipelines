#!/bin/bash

if [ $# -eq 0 ]
then
	echo "Usage: updatePipelineStatus.sh"
	echo "-q queueID"
	echo "-s status"
	exit
fi


if test $# -gt 0
	then
	while getopts :q:s: opt
	do
	case $opt in
	q)
		queueID=$OPTARG
		;;
	s)
			status=$OPTARG
		;;
	:)
		echo "Option -$OPTARG requires an argument."
		;;
	\?)
		echo "Invalid option: -$OPTARG"
	esac
	done
	shift $((OPTIND-1))
fi

username="root"
userpass="molSeq3127"
dbname="test"
tablename="pipelineStatus"

mysql -u $username -p$userpass -D $dbname -e "INSERT INTO $tablename
					(queueID,plStatus, timeUpdated) VALUES ('$queueID','$status', now());"

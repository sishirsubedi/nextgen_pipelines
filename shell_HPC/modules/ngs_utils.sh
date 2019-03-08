#!/bin/bash

log_info()
{
   log " [INFO] : $1"
}

log_error()
{
   log " [ERROR] : $1"
}

log()
{
 MESSAGE=$1
 TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
 SCRIPT=$( basename $0 )
 echo " [ $TIMESTAMP ] [ $SCRIPT ] $MESSAGE "
}

create_file()
{
  path=$1
  file_name=$2
	touch  ${path}${file_name}
	chmod 775  ${path}${file_name}
}

create_dir()
{
	new_dir=$1
	if [ ! -d $new_dir ] ; then
		mkdir $new_dir
	fi
	chmod 775 $new_dir
}

update_status()
{
user=$4
password=$5
database="ngs_$3"
insertstatement="INSERT INTO pipelineStatus (queueID, plStatus, timeUpdated) VALUES ('$1','$2',now());"
mysql --host="hhplabngsp01" --user="$user" --password="$password" --database="$database" --execute="$insertstatement"
}

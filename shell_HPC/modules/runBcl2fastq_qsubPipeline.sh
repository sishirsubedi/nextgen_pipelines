#!/bin/bash
#PBS -l mem=32gb,nodes=1
#PBS -l walltime=5:00:00
#PBS -q default


echo " from compute node: bcl2fastq is running for $DIR"

/usr/local/bin/bcl2fastq --no-lane-splitting --runfolder-dir $DIR --output-dir ${DIR}out2 -d 10 -p 10


wait $!


updatestatement="UPDATE pipelineStatusBcl2Fastq SET status=1 WHERE runID = $runID;"
mysql --user="$user" --password="$password" --database="$environment" --execute="$updatestatement"

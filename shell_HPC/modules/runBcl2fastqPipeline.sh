#!/bin/bash

run_bcl2fastq()
{
  USER=$1
	PASSWORD=$2
	DB=$3
	DB_HOST=$4
	RUNID=$5


	runFolder=$(ls -d /home/nextseq/*_"$RUNID"_*)


	qsub -k eo -F "-d$runFolder" /home/pipelines/bcl2fastq/runBcl2fastq_qsub.sh

}

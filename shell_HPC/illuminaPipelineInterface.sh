##############################################################################
#
# Houston Methodist Hospital
# Molecular Diagnostic
#
#Description:
#This script checks samples queued in instrument/assay specific txt file
# and calls appropriate interface script.
#Allocates appropriate PBS parameters.
##############################################################################

#!/bin/bash
#PBS -l nodes=1
#PBS -l walltime=2:00:00
#PBS -q default
#PBS -o ${PBS_JOBNAME}.out
#PBS -e ${PBS_JOBNAME}.err
#PBS -k eo

# ##############################################################################
# # functions
# ##############################################################################

display_usuage()
{
cat <<EOF >> /dev/stderr

 USAGE: $0

 OPTIONS:
 r - runID
 s - sampleName
 a - assay
 i - instrument
 e - environment
 q - queueID
 u - user
 p - password

EOF
}

parse_options()
{
    IMPORT=0

		while getopts "hz:r:s:a:i:e:q:u:p:" opt ; do
				case $opt in
					h)
						 display_usuage
						 exit 1
						 ;;
					z)
						IMPORT=$OPTARG
						;;
					r)
						RUNID=$OPTARG
						;;
					s)
					  SAMPLENAME=$OPTARG
						;;
					a)
						ASSAY=$OPTARG
						;;
					i)
				    INSTRUMENT=$OPTARG
					  ;;
				  e)
				    ENVIRONMENT=$OPTARG
					  ;;
				  q)
						QUEUEID=$OPTARG
						;;
					u)
						USER=$OPTARG
						;;
					p)
						PASSWORD=$OPTARG
						;;
					:)
						echo "Option -$OPTARG requires an argument."
						;;
					\?)
						echo "Invalid option: -$OPTARG"
			esac
	  done
		if [ $IMPORT -gt 0 ] ; then
				return 0
		fi
		return 1
}

load_modules()
{
      source /home/pipelines/ngs_${ENVIRONMENT}/shell/modules/ngs_utils.sh
}


create_rundate()
{
	##get runDate information##
	dateString=${RUNNAME%%_*}
  echo $dateString
	year=20${dateString:0:2}
	month=${dateString:2:2}
	day=${dateString:4:2}
	dateString=$year-$month-$day

	echo $dateString > ${HOME}runDate.txt

}

generate_ampliconFile(){

  HEME_EXCLUDED_DESIGN="/home/doc/ref/Heme/trusight-myeloid-amplicon-track.excluded.bed"

  /opt/samtools-1.4/samtools-1.4/samtools bedcov  $HEME_EXCLUDED_DESIGN  ${HOME}variantCaller/${SAMPLENAME}.sort.bam | awk ' {print $4,"\t",int($13/($8-$7))} ' > ${HOME}variantAnalysis/${SAMPLENAME}.samtools.coverageDepth

  AMPLICON_FILE=$(ls ${HOME}variantAnalysis/${SAMPLENAME}.samtools.coverageDepth)

}


prep_heme()
{


  HEME_EXCLUDED_DESIGN="/home/doc/ref/Heme/trusight-myeloid-amplicon-track.excluded.bed"
  #bedtools -u flag to write original A entry once if any overlaps found in B
  /opt/software/bedtools-2.17.0/bin/bedtools intersect -u -a $VARIANT_VCF -b $HEME_EXCLUDED_DESIGN > ${HOME}variantAnalysis/${SAMPLENAME}.filter.vcf


  HEME_EXCLUDED_AMPLICON="/home/doc/ref/Heme/excludedAmplicons.txt"
  # -v means "invert the match" in grep, in other words, return all non matching lines.
  grep -v -f $HEME_EXCLUDED_AMPLICON $AMPLICON_FILE > ${HOME}variantAnalysis/${SAMPLENAME}.amplicon.filter.txt


}


run_illuminaPipeline()
{
  bash ${HOME_SHELLDIR}illuminaPipeline.sh -d $HOME \
        -s $SAMPLENAME -e $ENVIRONMENT -q $QUEUEID -u $USER -p $PASSWORD
}
# ##############################################################################
# main
# ##############################################################################
main()
{
    parse_options $*

    if [ $? -eq 0 ]
    then
			  echo "Error in previous step. Aborting $0"
        exit 0
    fi

		############################################################################
		# initialize variables
		############################################################################

    load_modules


    RUNFOLDER=$(ls -d /home/$INSTRUMENT/*_"$RUNID"_*)
    RUNNAME=${RUNFOLDER##/home/$INSTRUMENT/}

    HOME="/home/environments/ngs_${ENVIRONMENT}/${INSTRUMENT}Analysis/${RUNNAME}/${SAMPLENAME}/"
		HOME_SHELLDIR="/home/pipelines/ngs_${ENVIRONMENT}/shell/"

    ##Aligning fastq files
    echo "Aligning fastq files"
    FASTQ_R1=$RUNFOLDER/out1/"$SAMPLENAME"*_R1_001.fastq.gz
    FASTQ_R2=${FASTQ_R1/_R1_/_R2_}      #replace "R1" with "R2"
    echo $FASTQ_R1
    echo $FASTQ_R2

    VARIANT_VCF=""
    AMPLICON_FILE=""

		create_dir ${HOME}variantCaller
		create_dir ${HOME}variantAnalysis

		exec >  >(tee -a ${HOME}process.log)
		exec 2> >(tee -a ${HOME}process.log >&2)

    show_pbsinfo

    log_info " Running illumina Pipeline Interface BY qsub for :
		assay : $ASSAY
		instrument : $INSTRUMENT
		runID : $RUNID
		sampleName : $SAMPLENAME
		environment : $ENVIRONMENT
		queueID : $QUEUEID "

    create_rundate

    if [ $INSTRUMENT == "miseq" ] ; then

      VARIANT_VCF=$(ls /home/${INSTRUMENT}/*_${RUNID}_*/Data/Intensities/BaseCalls/Alignment/${SAMPLENAME}_*.vcf)
      AMPLICON_FILE=$(ls /home/${INSTRUMENT}/*_${RUNID}_*/Data/Intensities/BaseCalls/Alignment/AmpliconCoverage_M1.tsv)

      echo $VARIANT_VCF
      echo $AMPLICON_FILE

    if [ $INSTRUMENT == "nextseq" ] ; then

	     bash ${HOME_SHELLDIR}VarScanPipelinePE.sh -p $FASTQ_R1 -q $FASTQ_R2 -o ${HOME}variantCaller -e $ENVIRONMENT


       VARIANT_VCF=$(ls ${HOME}variantCaller/${SAMPLENAME}.comb.vcf)
       generate_ampliconFile

       echo $VARIANT_VCF
       echo $AMPLICON_FILE

       prep_heme


    fi


    prep_heme


    #update_status "$QUEUEID" "Started" "$ENVIRONMENT" "$USER"  "$PASSWORD"

    run_illuminaPipeline



}


# ##############################################################################
# run main
# ##############################################################################


main $*

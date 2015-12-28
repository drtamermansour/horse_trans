#!/bin/sh

if [ $# -lt 4 ]
then
printf "\nUsage run_readGroupInfo_illumina.sh [samples list] [replicates list] [library Name] [script]\n"
exit 0
fi


sample_list="$1"
replicates_list="$2"
lib="$3"
script="$4"

module load SAMTools/0.1.19

while read bam; do
  sample_folder=$(dirname $bam | xargs basename)
  SM=${sample_folder#tophat_}                                                 ##smaple ID
  if [ -f $replicates_list ];then
    repLib_temp=$(grep $SM $replicates_list | awk '{ print $1 }')
    repLib=${repLib_temp%_R1_*.fastq.gz}
    if [ "$repLib" != "" ];then echo "found in replicates_list";
    else echo "Not found in replicates_list";fi
  else repLib=$SM;fi
  LB=$lib.$repLib                                                                     ##library ID
  PL="Illumina"                                                               ##platform (e.g. illumina, solid)
  ##read BAM 1st read, check the format.If typical, identify ID as "<instrument>:<run number>:<flowcell ID>:<lane>", and PU as the "<instrument>"
  header=$(samtools view $bam | head -n1 | awk '{ print $1}' | grep ':*:*:*:*:*:*')
  if [ "$header" != "" ]; then
    PU=$(echo ${header} | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)                ##platform unit-lane ID
  else # "make unique ID and PU using checksum"
    checksum=$(shasum $bam | awk '{ print $1 }')
    PU="UnChrID_"$checksum
  fi
  RGID=$PU.$SM

  output=${bam%.bam}_RG.bam
  cd $(dirname $bam)
  echo RGID $RGID LB $LB PL $PL PU $PU SM $SM
  qsub -v bam=$bam,output=$output,RGID=$RGID,LB=$LB,PL=$PL,PU=$PU,SM=$SM "${script}"
done < "$sample_list"

## info
# Description of the @RG items: The ID is usually a unique identifier for each lane. SM is the smaple ID. LB is the library identifier and PU refers to the sequencing machine.
# default fastq read name description is: @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<barcode sequence>
# https://www.broadinstitute.org/gatk/guide/best-practices?bpm=DNAseq#tutorials_mapdedup2799
# https://www.biostars.org/p/43897/
# https://www.biostars.org/p/47487/


#!/bin/sh

if [ $# -lt 3 ]
then
printf "\nUsage run_readGroupInfo_illumina.sh [samples list] [library Name] [script]\n"
exit 0
fi


sample_list="$1"
lib="$2"
script="$3"

while read bam; do
  ##read BAM 1st read, check the format.If typical, identify ID as "<instrument>:<run number>:<flowcell ID>:<lane>", and PU as the "<instrument>"
  header=$(samtools view $bam | head -n1 | awk '{ print $1}' | grep ':*:*:*:*:*:*')
  if [ "$header" != "" ]; then
    ID=$(echo ${header} | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)                ##lane ID
    PU=$(echo ${header} | cut -d ":" -f1)                                     ##platform unit (eg. run barcode)
  else # "make unique ID and PU using checksum"
    checksum=$(shasum $bam)
    ID="UnChrID_"$checksum
    PU="UnChrPU_"$checksum
  fi

  sample_folder=$(dirname $bam | xargs basename)
  SM=${sample_folder#tophat_}                                                 ##smaple ID
  LB=$lib                                                                     ##library ID
  PL="Illumina"                                                               ##platform (e.g. illumina, solid)
  output=${bam%.bam}_RG.bam
  cd $(dirname $bam)
  qsub -v bam=$bam,output=$output,ID=$ID,LB=$LB,PL=$PL,PU=$PU,SM=$SM "${script}"
done < "$sample_list"

## info
# Description of the @RG items: The ID is usually a unique identifier for each lane. SM is the smaple ID. LB is the library identifier and PU refers to the sequencing machine.
# default fastq read name description is: @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<barcode sequence>
# https://www.broadinstitute.org/gatk/guide/best-practices?bpm=DNAseq#tutorials_mapdedup2799
# https://www.biostars.org/p/43897/
# https://www.biostars.org/p/47487/


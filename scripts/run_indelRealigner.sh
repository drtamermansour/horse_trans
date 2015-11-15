#!/bin/sh

if [ $# -lt 6 ]
then
printf "\nUsage run_indelRealigner.sh [known indels] [indel intervals] [indexed reference fasta] [samples list] [the name of target BAM file] [script]\n"
exit 0
fi

indels="$1"
intervals="$2"
gatk_ref="$3"
sample_list="$4"
target_bam="$5"
script="$6"

while read f; do
  echo $f
  cd "$f"
  qsub -v indels="${indels}",intervals="${intervals}",gatk_ref="${gatk_ref}",sample=${target_bam} "${script}"
  cd ../
done < $sample_list
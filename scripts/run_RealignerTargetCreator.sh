#!/bin/sh

if [ $# -lt 5 ]
then
printf "\nUsage run_RealignerTargetCreator.sh [known indels] [indexed reference fasta] [samples list] [the name of target BAM file] [script]\n"
exit 0
fi

indels="$1"
gatk_ref="$2"
sample_list="$3"
target_bam="$4"
script="$5"

while read f; do
  echo $f
  cd "$f"
  qsub -v indels="${indels}",gatk_ref="${gatk_ref}",sample=${target_bam} "${script}"
  cd ../
done < $sample_list

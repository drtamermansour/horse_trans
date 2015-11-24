#!/bin/sh

if [ $# -lt 4 ]
then
printf "\nUsage run_applyRecalib.sh [indexed reference fasta] [samples list] [the name of target BAM file] [script]\n"
exit 0
fi

gatk_ref="$1"
sample_list="$2"
target_bam="$3"
script="$4"

while read f; do
  echo $f
  cd "$f"
  qsub -v gatk_ref="${gatk_ref}",sample="${target_bam}" "${script}"
  cd ../
done < $sample_list
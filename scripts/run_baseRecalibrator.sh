#!/bin/sh

if [ $# -lt 6 ]
then
printf "\nUsage run_baseRecalibrator.sh [known SNPs] [known indels] [indexed reference fasta] [samples list] [the name of target BAM file] [script]\n"
exit 0
fi

snps="$1"
indels="$2"
gatk_ref="$3"
sample_list="$4"
target_bam="$5"
script="$6"

while read f; do
  echo $f
  cd "$f"
  qsub -v snps="$snps",indels="${indels}",gatk_ref="${gatk_ref}",sample=${target_bam} "${script}"
  cd ../
done < $sample_list
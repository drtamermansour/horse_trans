#!/bin/sh

if [ $# -lt 5 ]
then
printf "\nUsage run_haplotypeCaller.sh [known SNPs] [indexed reference fasta] [samples list] [the name of target BAM file] [script]\n"
exit 0
fi

snps="$1"
gatk_ref="$2"
sample_list="$3"
target_bam="$4"
script="$5"

while read f; do
  echo $f
  cd "$f"
  qsub -v snps="$snps",gatk_ref="${gatk_ref}",sample=${target_bam} "${script}"
  cd ../
done < $sample_list
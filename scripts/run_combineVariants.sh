#!/bin/sh

if [ $# -lt 3 ]
then
printf "\nUsage run_haplotypeCaller.sh [indexed reference fasta] [samples list] [script]\n"
exit 0
fi

gatk_ref="$1"
sample_list="$2"
script="$3"

samples=""
while read f; do
  echo $f
  samples=" ${samples} --variant ${f}"
done < $sample_list

trim_samples=$(echo $samples | xargs | sed 's/\n//')
qsub -v gatk_ref="${gatk_ref}",samples="${trim_samples}" "${script}"

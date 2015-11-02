#!/bin/sh

if [ $# -lt 4 ]
then
printf "\nUsage run_indelRealigner.sh [known indels] [indexed reference fasta] [samples list] [script]\n"
exit 0
fi

indels="$1"
gatk_ref="$2"
sample_list="$3"
script="$4"

samples=""
while read f; do
  sample="$f"/$"split.bam"
  echo $sample;
  samples=" ${samples} -I ${sample}"
done < $sample_list

trim_samples=$(echo $samples | xargs)
qsub -v indels="${indels}",gatk_ref="${gatk_ref}",samples="${trim_samples}" "${script}"

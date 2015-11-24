
#!/bin/sh

if [ $# -lt 4 ]
then
printf "\nUsage run_haplotypeCaller.sh [indexed reference fasta] [samples list] [the name of target BAM file] [script]\n"
exit 0
fi

gatk_ref="$1"
sample_list="$2"
target_bam="$3"
script="$4"

samples=""
while read f; do
  echo $f
  sample="$f"/"$target_bam"
  samples=" ${samples} -I ${sample}"
done < $sample_list

trim_samples=$(echo $samples | xargs | sed 's/\n//')
qsub -v gatk_ref="${gatk_ref}",samples="${trim_samples}" "${script}"

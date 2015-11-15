#!/bin/sh

if [ $# -lt 5 ]
then
printf "\nUsage run_baseRecalibrator.sh [list of known variants] [indexed reference fasta] [samples list] [the name of target BAM file] [script]\n"
exit 0
fi

known_var="$1"
gatk_ref="$2"
sample_list="$3"
target_bam="$4"
script="$5"


variants=""
while read var; do
  variants=" ${variants} -knownSites ${var}"
done < $known_var
trim_variants=$(echo $variants | xargs)
echo $trim_variants


#while read f; do
#  echo $f
#  cd "$f"
#  qsub -v variants="${trim_variants}",gatk_ref="${gatk_ref}",sample="${target_bam}" "${script}"
#  cd ../
#done < $sample_list

samples=""
while read f; do
  sample="$f"/${target_bam}
  samples=" ${samples} -I ${sample}"
done < $sample_list
trim_samples=$(echo $samples | xargs)
echo $samples;

qsub -v variants="${trim_variants}",gatk_ref="${gatk_ref}",samples="${trim_samples}" "${script}"



















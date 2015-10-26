#!/bin/sh

if [ $# -lt 2 ]
then
    printf "\nUsage run_cufflinks_noRef.sh [samples list] [script]\n"
    exit 0
fi

sample_list="$1"
script="$2"


while read f; do
  cd "$f"
  f2=$(basename $f)
  label=$"Sample_"${f2#tophat_}
  sample=$"accepted_hits_RG.bam"
  if [ ! -f "$sample" ];then
    sample=$"accepted_hits.bam"; fi
  echo $label $sample;
  qsub -v label=${label},sample=${sample} "${script}"
  cd ..
done < $sample_list

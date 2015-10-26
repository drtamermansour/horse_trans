#!/bin/sh

if [ $# -lt 2 ]
then
printf "\nUsage run_BuildBamIndex.sh [samples list] [script]\n"
exit 0
fi

sample_list="$1"
script="$2"


while read f; do
  cd "$f"
  #f2=$(basename $f)
  #label=$"Sample_"${f2#tophat_}
  sample=$"reordered.bam"
  echo $sample;
  qsub -v sample=${sample} "${script}"
  cd ..
done < $sample_list

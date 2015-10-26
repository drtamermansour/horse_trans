#!/bin/sh

if [ $# -lt 2 ]
then
printf "\nUsage run_reorderBAM.sh [samples list] [script]\n"
exit 0
fi

gatk_ref="$1"
sample_list="$2"
script="$3"


while read f; do
cd "$f"
#f2=$(basename $f)
#label=$"Sample_"${f2#tophat_}
sample=$"dedup_reads.bam"
echo $sample;
qsub -v gatk_ref=${gatk_ref},sample=${sample} "${script}"
cd ..
done < $sample_list

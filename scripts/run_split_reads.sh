#!/bin/sh

if [ $# -lt 2 ]
then
printf "\nUsage run_split_reads.sh [work_dir] [script]\n"
exit 0
fi


sample_list="$1"
script="$2"

while read f; do echo $f;
  base=$(basename ${f%_R1_*.pe.se.fq})
  fkeep="$base"_R_*.pe.fq.keep
  f1="$base"_R1_001.pe.fq
  f2="$base"_R2_001.pe.fq
  qsub -v fkeep=$fkeep,f1=$f1,f2=$f2 $script
done < $sample_list

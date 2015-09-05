#!/bin/sh

if [ $# -lt 3 ]; then
  printf "\nUsage run_cufflinks.sh [samples list] [library type] [script]\n"
  exit 0
fi

sample_list="$1"
lib="$2"
script="$3"


while read f; do
  basef=$(basename $f)
  sample=${basef#tophat_}
  pe=$sample"_R_001.pe.fq"
  se=$sample"_R_001.se.fq"
  qsub -v input_dir=$f,sample=$sample,output_pe=$pe,output_se=$se,lib_type=$lib $script
done < $sample_list




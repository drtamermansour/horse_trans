#!/bin/sh

if [ $# -lt 3 ]
then
    printf "\nUsage run_cufflinks_wRef.sh [samples list] [Genes_GTF_file] [script]\n"
    exit 0
fi

sample_list="$1"
Genes_GTF_file="$2"
script="$3"


while read f; do
  cd "$f"
  f2=$(basename $f)
  label=$"Sample_"${f2#tophat_}
  sample=$"accepted_hits_RG.bam"
  if [ ! -f "$sample" ];then
    sample=$"accepted_hits.bam"; fi
  echo $label $sample;
  qsub -v Genes_GTF_file=${Genes_GTF_file},label=${label},sample=${sample} "${script}"
  cd ..
done < $sample_list

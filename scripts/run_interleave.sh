#!/bin/sh

if [ $# -lt 2 ]
then
printf "\nUsage run_tophat.sh [samples list] [script]\n"
exit 0
fi

sample_list="$1"
script="$2"


while read f; do echo $f;
  f1=$(echo $f | sed 's/.pe.se.fq/.pe.fq/');
  f2=$(echo $f1 | sed 's/_R1_/_R2_/');
  newf=$(echo $f1 | sed 's/_R1_/_R_/');
  qsub -v f1=$f1,f2=$f2,output=$newf $script
done < $sample_list



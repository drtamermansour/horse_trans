#!/bin/sh

if [ $# -lt 4 ]
then
  printf "\nUsage run_cuffcompare.sh [sample] [label] [log] [script]\n"
exit 0
fi

sample="$1"
label="$2"
log="$3"
script="$4"

qsub -v sample="${sample}",label="${label}",log="${log}" "${script}"


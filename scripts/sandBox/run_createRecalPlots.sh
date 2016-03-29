#!/bin/sh

if [ $# -lt 2 ]
then
printf "\nUsage run_createRecalPlots.sh [indexed reference fasta] [script]\n"
exit 0
fi

gatk_ref="$1"
script="$2"

qsub -v gatk_ref="${gatk_ref}" "${script}"

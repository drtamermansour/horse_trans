#!/bin/sh

if [ $# -lt 3 ]
then
printf "\nUsage run_cufflinks_wRef_single.sh [A sample] [sample label] [script]\n"
exit 0
fi

sample="$1"
label="$2"
script="$3"

qsub -v label=${label},sample=${sample} "${script}"

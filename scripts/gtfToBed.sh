#!/bin/sh

if [ $# -lt 3 ]
then
printf "\nUsage gtfToBed.sh [GTF file] [script path]\n"
exit 0
fi

targetAss="$1"
script_path="$2"


module load ucscUtils/262   # UCSC_kent_commands
if [ -f $targetAss ]; then
  filename=${targetAss%.gtf}
  gtfToGenePred $targetAss ${filename}.gpred
  cat ${filename}.gpred | $script_path/genePredToBed > ${filename}.bed
else
  printf "\nNo GTF found"; exit 0;
fi
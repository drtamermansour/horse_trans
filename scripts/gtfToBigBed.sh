#!/bin/sh

if [ $# -lt 3 ]
then
printf "\nUsage run_cufflinks.sh [genome chrom sizes] [identifier] [script path]\n"
exit 0
fi

chrom_sizes="$1"
identifier="$2"
script_path="$3"


module load ucscUtils/262   # UCSC_kent_commands
n=$(ls *.gtf | wc -l)
if [ $n -eq 1 ]; then
  targetAss=*.gtf
  filename=$(echo $identifier | sed 's/\//_/g' | sed 's/_output//g')
  gtfToGenePred $targetAss ${filename}.gpred
  cat ${filename}.gpred | $script_path/genePredToBed > ${filename}.bed
  sort -k1,1 -k2,2n ${filename}.bed > ${filename}_sorted.bed
  bedToBigBed ${filename}_sorted.bed $chrom_sizes ${filename}.BigBed
  #module load BEDOPS/2.4.5
  #gtf2bed < ${filename}.gtf > ${filename}.bedops
elif [ $n -eq 0 ]; then
  printf "\nNo GTF found"; exit 0;
elif [ $n -gt 1 ]; then
  printf "Many GTF files. Can not continue."; exit 0;
fi
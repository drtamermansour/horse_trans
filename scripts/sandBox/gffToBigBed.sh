#!/bin/sh

if [ $# -lt 3 ]
then
printf "\nUsage gffToBigBed.sh [genome chrom sizes] [identifier] [script path]\n"
exit 0
fi

chrom_sizes="$1"
identifier="$2"
script_path="$3"


#module load ucscUtils/262   # UCSC_kent_commands
n=$(ls *.gffp | wc -l)
if [ $n -eq 1 ]; then
  targetAss=*.gffp
  filename=$(echo $identifier | sed 's/\//_/g' | sed 's/_output//g')
  $script_path/UCSC_kent_commands/gff3ToGenePred $targetAss ${filename}.gpred
  cat ${filename}.gpred | $script_path/genePredToBed > ${filename}.bed
  sort -k1,1 -k2,2n ${filename}.bed > ${filename}_sorted.bed
  $script_path/UCSC_kent_commands/bedToBigBed ${filename}_sorted.bed $chrom_sizes ${filename}.BigBed
elif [ $n -eq 0 ]; then
  printf "\nNo The target GFF was not found"; exit 0;
elif [ $n -gt 1 ]; then
  printf "Many genome.gff3 files. Can not continue."; exit 0;
fi

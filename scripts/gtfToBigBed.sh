#!/bin/sh

if [ $# -lt 3 ]
then
printf "\nUsage gtfToBigBed.sh [GTF file] [genome chrom sizes] [script path]\n"
exit 0
fi

targetAss="$1"
chrom_sizes="$2"
script_path="$3"


#module load ucscUtils/262   # UCSC_kent_commands
if [ -f $targetAss ]; then
  filename=${targetAss%.gtf}
  $script_path/UCSC_kent_commands/gtfToGenePred $targetAss ${filename}.gpred
  cat ${filename}.gpred | $script_path/genePredToBed > ${filename}.bed
  sort -k1,1 -k2,2n ${filename}.bed > ${filename}_sorted.bed
  $script_path/UCSC_kent_commands/bedToBigBed ${filename}_sorted.bed $chrom_sizes ${filename}.BigBed
else
  printf "\nNo GTF found"; exit 0;
fi

#!/bin/sh

if [ $# -lt 2 ]
then
printf "\nUsage bedToBigBed.sh [BED file] [genome chrom sizes]\n"
exit 0
fi

targetAss="$1"
chrom_sizes="$2"
script_path=$(dirname "${BASH_SOURCE[0]}") 


#module load ucscUtils/262   # UCSC_kent_commands
if [ -f $targetAss ]; then
  filename=${targetAss%.bed}
  tail -n +2 $targetAss | sort -k1,1 -k2,2n > ${filename}_sorted.bed
  $script_path/UCSC_kent_commands/bedToBigBed ${filename}_sorted.bed $chrom_sizes ${filename}.BigBed
else
  printf "\nNo The target GFF was not found"; exit 0;
fi

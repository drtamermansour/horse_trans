#!/bin/sh

if [ $# -lt 2 ]
then
  printf "\nUsage run_tophat.sh [horse index] [transcriptome]\n"
  exit 0
fi

index="$1"
transcriptome="$2"


module load GNU/4.8.2
# The following have been reloaded with a version change:
# 1) Boost/1.47.0 => Boost/1.55.0     2) GNU/4.4.5 => GNU/4.8.2     3) OpenMPI/1.4.3 => OpenMPI/1.6.5     4) R/2.15.1 => R/3.1.0
module load salmon/0.5.0
# The following have been reloaded with a version change:
# 1) CMake/2.8.5 => CMake/3.1.0
salmon index --index "$index" --transcripts "$transcriptome" --type quasi


#!/bin/sh

if [ $# -lt 2 ]
then
  printf "\nUsage run_salmonIndex.sh [index name] [transcriptome]\n"
  exit 0
fi

index="$1"
transcriptome="$2"
script="$3"

qsub -v index=$index,transcriptome=$transcriptome $script


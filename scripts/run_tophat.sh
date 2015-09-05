#!/bin/sh

if [ $# -lt 6 ]
then
    printf "\nUsage run_tophat.sh [samples list] [library type] [first strand] [Bowtie2_genome_index_base] [transcriptome_index] [script_path]\n"
    exit 0
fi

sample_list="$1"
lib="$2"
strand="$3"
Bowtie2_genome_index_base="$4"
transcriptome_index="$5"
script_path="$6"

if [ $lib = $"PE" ]; then
  echo "Running PE mode";
  while read f; do
    base=$(basename "$f")
    base2=${base%_R1_*.pe.se.fq}
    output_dir=$"tophat_"$base2
    temp=$(echo $f | sed 's/_R1_/_R2_/')
    f2=$(echo $temp | sed 's/.pe.se.fq/.pe.fq/')
    echo $output_dir
    mkdir ${output_dir}
    qsub -v output_dir="${output_dir}",\
Bowtie2_genome_index_base="${Bowtie2_genome_index_base}",\
transcriptome_index="${transcriptome_index}",\
lib_type="${strand}",\
f1="$f",f2="$f2" $script_path/tophat_PE.sh
  done < $sample_list
else echo "No PE run"; fi

if [ $lib = $"SE" ]; then
  echo "Running SE mode";
  while read f; do
    base=$(basename "$f")
    base2=${base%_SR_*.se.fq}
    output_dir=$"tophat_"$base2
    echo $output_dir
    mkdir ${output_dir}
    qsub -v output_dir="${output_dir}",\
Bowtie2_genome_index_base="${Bowtie2_genome_index_base}",\
transcriptome_index="${transcriptome_index}",\
lib_type="${strand}",\
f1="$f" $script_path/tophat_SE.sh
  done < $sample_list
else echo "No SE run"; fi



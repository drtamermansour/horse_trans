#!/bin/sh

if [ $# -lt 3 ]
then
    printf "\nUsage run_adapter_trimmer.sh [Sample List] [Library type] [platform]\n"
    exit 0
fi

sample_list="$1"
lib="$2"
platform="$3"
script_path=$(dirname "${BASH_SOURCE[0]}")  

if [ $lib = $"PE" ]; then
  while read f; do
    echo $f
    input_one="$f"
    input_two=$(echo "$f" | sed s/_R1_/_R2_/)
    temp=$(basename "$f")
    temp2=${temp%.fastq.gz}
    output_pe1=$temp2".pe.fq"
    output_pe2=$(echo "$output_pe1" | sed s/_R1_/_R2_/)
    output_se1=$temp2".se.fq"
    output_se2=$(echo "$output_se1" | sed s/_R1_/_R2_/)
    if [ $platform == $"HPC" ];then
      qsub -v R1_INPUT="$input_one",R2_INPUT="$input_two",output_pe1="$output_pe1",output_pe2="$output_pe2",output_se1="$output_se1",output_se2="$output_se2" $script_path/adapter_trimmer_PE.sh
    elif [ $platform == $"AMZ" ];then
      bash $script_path/adapter_trimmer_PE_AMZ.sh "$input_one" "$input_two" "$output_pe1" "$output_pe2" "$output_se1" "$output_se2" 
    else echo "Platform should be either HPC or AMZ";fi
  done < $sample_list

elif [ $lib = $"SE" ]; then
  while read f; do
    echo $f
    temp=$(basename "$f")
    temp2=${temp%.fastq.gz}
    output=$temp2".se.fq"
    if [ $platform == $"HPC" ];then 
      qsub -v INPUT="$f",output="$output" $script_path/adapter_trimmer_SE.sh
    elif [ $platform == $"AMZ" ];then 
    ## prep adapter_trimmer_SE_AMZ.sh
    else echo "Platform should be either HPC or AMZ";fi
  done < $sample_list

else echo $"Inappropriate format of target directory"; fi;







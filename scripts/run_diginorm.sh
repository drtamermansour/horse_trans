#!/bin/sh

if [ $# -lt 4 ]
then
printf "\nUsage run_diginorm.sh [work_dir] [kmer] [cutoff] [script_path]\n"
exit 0
fi


work_dir="$1"
kmer="$2"
cutoff="$3"
script_path="$4"

seq=$(basename $work_dir | cut -d"_" -f1)

kmer_table=$"C$cutoff"k"$kmer".kt

PEsize=$(du $work_dir/trimmed_RNA_reads/*_R_*.pe.fq | awk '{ total += $1 }; END { rounded = sprintf("%.0f", total/1000000); print rounded }')
SEsize=$(du $work_dir/trimmed_RNA_reads/*_*R_*.se.fq | awk '{ total += $1 }; END { rounded = sprintf("%.0f", total/1000000); print rounded }')
total_size=$(echo "$PEsize + $SEsize" | bc)

if [ $total_size -lt 210 ]; then
  threadRAM=$(echo "$total_size/10 + 11" | bc)
  walltime=$(echo "$total_size/10 + 4" | bc)
elif [ $total_size -lt 1650 ]; then
  threadRAM=31; walltime=$(echo "($total_size/10 + 3)/24" | bc)":00";
else threadRAM=31; walltime=$"07:00"; fi

RAM=$(echo "$threadRAM*4 + 4" | bc)"Gb"

if [ $seq = $"PE" ]; then
  echo "Running diginorm in paired mode";
  cat ${work_dir}/trimmed_RNA_reads/*_*R_*.se.fq > ${work_dir}/trimmed_RNA_reads/allsingletons.fq
  samples=()
  for f in ${work_dir}/trimmed_RNA_reads/*_R_*.pe.fq; do
    samples+=($f); done;
  qsub -l walltime="$walltime":00:00,nodes=1:ppn=1,mem="$RAM" \
-v kmer=$kmer,cutoff=$cutoff,x=$threadRAM"e9",kmer_table=$kmer_table,\
singletons="${work_dir}/trimmed_RNA_reads/allsingletons.fq",\
samples="${samples[*]}" $script_path/diginorm_pe.sh
else echo "No paired end run"; fi

if [ $seq = $"SE" ]; then
  echo "Running diginorm in single end mode";
  cat ${work_dir}/trimmed_RNA_reads/*_*R_*.se.fq > ${work_dir}/trimmed_RNA_reads/allsingletons.fq
  qsub -l walltime="$walltime":00:00,nodes=1:ppn=1,mem="$RAM" \
-v kmer=$kmer,cutoff=$cutoff,x=$threadRAM"e9",kmer_table=$kmer_table,\
singletons="${work_dir}/trimmed_RNA_reads/allsingletons.fq" $script_path/diginorm_se.sh
else echo "No single end run"; fi



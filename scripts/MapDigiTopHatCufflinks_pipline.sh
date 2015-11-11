#!/bin/sh
myRoot=$"/mnt/ls12/Tamer"
source $myRoot/config.txt
###########################################################################################
## pipeline_mapped_diginormAllsamples_Tophat2.refGTFguided_Cufflinks.refGTFguided.Cuffmerge
## convert the tophat output BAM files into Fastq files
while read work_dir; do
  echo $work_dir
  mkdir -p $work_dir/mapped_trimmed
  sample_list=$work_dir/tophat_output/sample_list.txt
  lib=$(basename $work_dir | cut -d"_" -f 1)
  bash ${script_path}/run_BamToFastq.sh "$sample_list" "$lib" "$work_dir/mapped_trimmed" "$script_path/restore_mapped_trimmed.sh"
done < $horse_trans/working_list_NoPBMCs_NoCereb.txt

## Check for successful jobs (requires restore_mapped.e & restore_mapped.o)
while read work_dir; do
  cd $work_dir/mapped_trimmed
  sample_list=$work_dir/mapped_trimmed/failedSamples.txt           ## define the path of empty file
  bash $script_path/check_BamToFastq.sh "$sample_list"
  x=$(cat $sample_list | wc -l)
  if [ $x -ne 0 ]; then
    lib=$(basename $work_dir | cut -d"_" -f 1)
    echo "Failed jobs in: "$work_dir
    bash ${script_path}/run_BamToFastq.sh "$sample_list" "$lib" "$script_path/restore_mapped_trimmed.sh"
fi; done < $horse_trans/working_list_NoPBMCs_NoCereb.txt

## run digital normalization of lab specific tissues (need to be updated to use sample list and check for success)
kmer=20
cutoff=50
while read work_dir; do
  echo $work_dir
  mkdir -p $work_dir/normalizied_mappedRNA_reads
  cd $work_dir/normalizied_mappedRNA_reads
  lib=$(basename $work_dir | cut -d"_" -f 1)
  bash ${script_path}/run_diginorm.sh "$lib" "$work_dir/mapped_trimmed" "$kmer" "$cutoff" "$script_path"
done < $horse_trans/working_list_NoPBMCs_NoCereb.txt

## Check for successful diginorm and trouble shooting the failed jobs (requires T_Trim.e)
sample_list=$prepData/failed_diginorm.txt                   ## define the path of empty file
> $sample_list
while read work_dir; do
  cd $work_dir/normalizied_mappedRNA_reads
  bash $script_path/check_diginorm.sh "$sample_list"
done < $horse_trans/working_list_NoPBMCs_NoCereb.txt
x=$(cat $sample_list | wc -l)
if [ $x -ne 0 ]; then
    echo "Failed jobs in: "
    cat $sample_list
    while read work_dir; do
      cd $work_dir/normalizied_mappedRNA_reads
      lib=$(basename $work_dir | cut -d"_" -f 1)
      bash ${script_path}/run_diginorm.sh "$lib" "$work_dir/mapped_trimmed" "$kmer" "$cutoff" "$script_path"
    done < $sample_list
fi

## split the interleaved reads
while read work_dir; do
  echo $work_dir
  cd $work_dir/normalizied_mappedRNA_reads
  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
  if [ "$lib" = $"PE" ]; then
    #sample_list=$work_dir/trimmed_RNA_reads/sample_list.txt
    sample_list=$work_dir/tophat_output/sample_list.txt
    bash ${script_path}/run_split_reads2.sh "$sample_list" $script_path/split_reads.sh
fi; done < $horse_trans/working_list_NoPBMCs.txt

## merge singletons and change the file names to fir the tophat script
while read work_dir; do
  echo $work_dir
  cd $work_dir/normalizied_mappedRNA_reads
  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
  if [ "$lib" = $"PE" ]; then
    singletones=1
    for f in *_R1_001.pe.fq; do
      base=${f%_R1_001.pe.fq}
      if [ $singletones -eq 1 ]; then cat $f allsingletons.fq.keep > "$base"_R1_001.pe.se.fq; singletones=0;
      else mv $f "$base"_R1_001.pe.se.fq; fi; done
  elif [ "$lib" = $"SE" ]; then
    mv allsingletons.fq.keep allsingletons_SR_002.se.fq
fi; done < $horse_trans/working_list_NoPBMCs.txt

## define the list samples.
## This is where you can edit the output list file(s) to restrict the processing for certain target(s)
while read work_dir; do if [ -d $work_dir/normalizied_RNA_reads ]; then
rm -f $work_dir/normalizied_RNA_reads/sample_list.txt
for f in $work_dir/normalizied_RNA_reads/{*_R1_*.pe.se.fq,*_SR_*.se.fq}; do if [ -f $f ]; then
echo $f >> $work_dir/normalizied_RNA_reads/sample_list.txt; fi; done;
fi; done < $horse_trans/working_list_NoPBMCs.txt

## run Tophat on each sample
while read work_dir; do
echo $work_dir
mkdir -p $work_dir/digi_tophat_output
cd $work_dir/digi_tophat_output
lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
strand=$(basename $work_dir | cut -d"_" -f 3 | sed 's/\./-/')   ## fr-unstranded, fr-firststrand or fr-secondstrand
sample_list=$work_dir/normalizied_RNA_reads/sample_list.txt
bash ${script_path}/run_tophat.sh "$sample_list" "$lib" "$strand" "$Bowtie2_genome_index_base" "$transcriptome_index" "$script_path"
done < $horse_trans/working_list_NoPBMCs.txt










#################################################################################
## Further analysis for unmapped reads
while read work_dir; do
  echo $work_dir
  mkdir -p $work_dir/unmapped_trimmed
  cd $work_dir/unmapped_trimmed
  sample_list=$work_dir/tophat_output/sample_list.txt
  lib=$(basename $work_dir | cut -d"_" -f 1)
  bash ${script_path}/run_BamToFastq.sh "$sample_list" "$lib" "$script_path/restore_unmapped_trimmed.sh"
done < $horse_trans/working_list_NoPBMCs_NoCereb_NoSkin.txt


#!/bin/sh
myRoot=$"/mnt/ls12/Tamer"
source $myRoot/config.txt
###########################################################################################
## pipeline_diginormAllsamples_Tophat2.refGTFguided_Cufflinks.refGTFguided.Cuffmerge
cd ~/khmer
git checkout horseTrans         ## the branch was created on 08/24/2015 for reproducibility
## interleave PE files
while read work_dir; do
  echo $work_dir
  cd $work_dir/trimmed_RNA_reads
  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
  if [ "$lib" = $"PE" ]; then
    sample_list=$work_dir/trimmed_RNA_reads/sample_list.txt
    bash ${script_path}/run_interleave.sh "$sample_list" $script_path/interleave.sh
fi; done < $horse_trans/working_list_NoPBMCs.txt

## Check for successful interleave runs and trouble shooting the failed jobs (requires interleave.e)
while read work_dir; do
  cd $work_dir/trimmed_RNA_reads
  sample_list=$work_dir/trimmed_RNA_reads/interleave_failedSamples.txt           ## define the path of empty file
  bash $script_path/check_interleave.sh "$sample_list"
  x=$(cat $sample_list | wc -l)
  if [ $x -ne 0 ]; then
    echo "Failed interleave jobs in: "$work_dir
    cat $sample_list
    bash ${script_path}/run_interleave.sh "$sample_list" $script_path/interleave.sh;
fi
done < $horse_trans/working_list_NoPBMCs.txt

## run digital normalization of lab specific tissues (need to be updated to use sample list and check for success)
kmer=20
cutoff=50
while read work_dir; do
  echo $work_dir
  mkdir -p $work_dir/normalizied_RNA_reads
  cd $work_dir/normalizied_RNA_reads
  lib=$(basename $work_dir | cut -d"_" -f 1)
  bash ${script_path}/run_diginorm.sh "$lib" "$work_dir/trimmed_RNA_reads" "$kmer" "$cutoff" "$script_path"
done < $horse_trans/working_list_NoPBMCs.txt

## Check for successful diginorm and trouble shooting the failed jobs (requires T_Trim.e)
sample_list=$prepData/failed_diginorm.txt                   ## define the path of empty file
> $sample_list
while read work_dir; do
  cd $work_dir/normalizied_RNA_reads
  bash $script_path/check_diginorm.sh "$sample_list"
done < $horse_trans/working_list_NoPBMCs.txt
x=$(cat $sample_list | wc -l)
if [ $x -ne 0 ]; then
  echo "Failed jobs in: "
  cat $sample_list
  while read work_dir; do
    cd $work_dir/normalizied_RNA_reads
    lib=$(basename $work_dir | cut -d"_" -f 1)
    bash ${script_path}/run_diginorm.sh "$lib" "$work_dir/trimmed_RNA_reads" "$kmer" "$cutoff" "$script_path"
  done < $sample_list
fi

## split the interleaved reads
while read work_dir; do
  echo $work_dir
  cd $work_dir/normalizied_RNA_reads
  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
  if [ "$lib" = $"PE" ]; then
    sample_list=$work_dir/trimmed_RNA_reads/sample_list.txt
    bash ${script_path}/run_split_reads.sh "$sample_list" $script_path/split_reads.sh
fi; done < $horse_trans/working_list_NoPBMCs.txt

## merge singletons and change the file names to fit the tophat script
while read work_dir; do
  echo $work_dir
  cd $work_dir/normalizied_RNA_reads
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

## merge the files in tophat compatible format
#while read work_dir; do
#  echo $work_dir
#  cd $work_dir/normalizied_RNA_reads
#  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
#  if [ "$lib" = $"PE" ]; then
#    cat *_R1_001.pe.fq allsingletons.fq.keep > allsamples_R1_002.pe.se.fq
#    cat *_R2_001.pe.fq > allsamples_R2_002.pe.fq
#  elif [ "$lib" = $"SE" ]; then
#    mv allsingletons.fq.keep allsingletons_SR_002.se.fq
#fi; done < $horse_trans/working_list_NoPBMCs.txt

## define the list samples.
## This is where you can edit the output list file(s) to restrict the processing for certain target(s)
while read work_dir; do if [ -d $work_dir/normalizied_RNA_reads ]; then
  rm -f $work_dir/normalizied_RNA_reads/sample_list.txt
  for f in $work_dir/normalizied_RNA_reads/{*_R1_*.pe.se.fq,*_SR_*.se.fq}; do if [ -f $f ]; then
    echo $f >> $work_dir/normalizied_RNA_reads/sample_list.txt; fi; done;
fi; done < $horse_trans/working_list_NoPBMCs.txt

## run Tophat on each sample independantly
## we meight need to run tophat on one sample after another feeding each run by the junctions of the previous run
while read work_dir; do
  echo $work_dir
  mkdir -p $work_dir/digi_tophat_output
  cd $work_dir/digi_tophat_output
  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
  strand=$(basename $work_dir | cut -d"_" -f 3 | sed 's/\./-/')   ## fr-unstranded, fr-firststrand or fr-secondstrand
  sample_list=$work_dir/normalizied_RNA_reads/sample_list.txt
  bash ${script_path}/run_tophat.sh "$sample_list" "$lib" "$strand" "$Bowtie2_genome_index_base" "$transcriptome_index" "$script_path"
done < $horse_trans/working_list_NoPBMCs.txt

## Check for successful tophat runs and trouble shooting the failed tophat jobs (require tophat-[SP]E.e & .o)
while read work_dir; do
  cd $work_dir/digi_tophat_output
  sample_list=$work_dir/normalizied_RNA_reads/tophat_failedSamples.txt           ## define the path of empty file
  bash $script_path/check_tophat.sh "$sample_list"
  x=$(cat $sample_list | wc -l)
  if [ $x -ne 0 ]; then
    lib=$(basename $work_dir | cut -d"_" -f 1)
    strand=$(basename $work_dir | cut -d"_" -f 3 | sed 's/\./-/')
    echo "Failed tophat jobs in: "$work_dir
    bash ${script_path}/run_tophat.sh "$sample_list" "$lib" "$strand" "$Bowtie2_genome_index_base" "$transcriptome_index" "$script_path"
fi; done < $horse_trans/working_list_NoPBMCs.txt

##################
## create summary for tophat run
while read work_dir; do
  for f in $work_dir/digi_tophat_output/tophat_*; do
    echo ${f} >> $work_dir/digi_tophat_output/allsample_summary.txt
    cd ${f}
    grep "overall read mapping rate" align_summary.txt >> ../allsample_summary.txt
    grep "concordant pair alignment rate" align_summary.txt >> ../allsample_summary.txt
  done
done < $horse_trans/working_list_NoPBMCs.txt
##################
## define the list samples.
## This is where you can edit the output list file(s) to restrict the processing for certain target(s)
while read work_dir; do if [ -d $work_dir/digi_tophat_output ]; then
  rm -f $work_dir/digi_tophat_output/sample_list.txt
  for f in $work_dir/digi_tophat_output/tophat_*; do if [ -d $f ]; then
    echo $f >> $work_dir/digi_tophat_output/sample_list.txt; fi; done;
fi; done < $horse_trans/working_list_NoPBMCs.txt

## Merge BAM files
module load SAMTools/0.1.19
while read work_dir; do
  echo $work_dir
  cd $work_dir/digi_tophat_output
  samples=()
  while read sample; do samples+=($sample/accepted_hits.bam); done < sample_list.txt
  len=${#samples[@]}
  if [ $len -gt 1 ]; then samtools merge merged.bam ${samples[*]};
  elif [ $len -eq 1 ]; then cp ${samples[0]} merged.bam;
  else echo "can find bam files"; fi
done < $horse_trans/working_list_NoPBMCs.txt

### Run Cufflinks: output transcripts.gtf in the same tophat_sample folder
while read work_dir; do
  echo $work_dir
  cd $work_dir/digi_tophat_output
  sample=$work_dir/digi_tophat_output/merged.bam
  label=$(basename $work_dir)
  bash ${script_path}/run_cufflinks2.sh "$sample" "$refGTF_file" "$label" "$script_path/cufflinks2.sh";
done < $horse_trans/working_list_NoPBMCs.txt

## Check for successful Cufflinks runs and trouble shooting the failed Cufflinks jobs (requires cufflinks.e)
while read work_dir; do if [ -d $work_dir/digi_tophat_output ]; then
  echo $work_dir
  cd $work_dir/digi_tophat_output
  sample_list=$work_dir/digi_tophat_output/failedSamples.txt           ## define the path of empty file
  bash $script_path/check_cufflinks2.sh "$sample_list"
  x=$(cat $sample_list | wc -l)
  if [ $x -ne 0 ]; then
    echo "Failed Cufflinks jobs in: "$work_dir
    sample=$(cat $sample_list)
    label=$(basename $work_dir)
    bash ${script_path}/run_cufflinks2.sh "$sample" "$refGTF_file" "$label" "$script_path/cufflinks2.sh";
fi; fi; done < $horse_trans/working_list_NoPBMCs.txt

# For every target tissue, diginorm the lab-spefific BAMs $tissue_merge/digimerge/TISSUE.NAME/withORwithoutRefGuidence.gtf
# 1. diginorm: start with tissues with longer reads (& higher mapping effeciency)
# 2. tophat libraries separately
# 3. merge BAM files
# 4. Cufflinks
# 5. ln to the tissue_merge path
tissue_Digimerge=$tissue_merge/digimerge
mkdir -p $tissue_Digimerge
cd $tissue_Digimerge
while read tissue_dir; do
  tissue=$(basename $tissue_dir)
  mkdir -p $tissue
  ln XXXXXXXXXXXXXXXXXXXXXX
done < $horse_trans/multi_lib_tissues.txt
##################
## multi-tissue diginorm into one total assembly
# save the assembly in $horse_trans/total_merge/cuffmerge/withORwithoutRefGuidence.gtf

###################
## create list of assemblies from each library
## This is where you can edit the list to restrict the processing for certain target(s)
rm -f $prepData/digi_merged_assemblies.txt
rm -f $horse_trans/rawdigi_TopCuff_assemblies.txt
while read work_dir; do if [ -d $work_dir/digi_tophat_output ]; then
  dir=$work_dir/digi_tophat_output
  mkdir $dir/withoutGTF
  ln $dir/transcripts.gtf $dir/withoutGTF/.
  echo ${dir#$prepData/}/withoutGTF >> $prepData/digi_merged_assemblies.txt;
  echo ${dir}/withoutGTF >> $horse_trans/rawdigi_TopCuff_assemblies.txt; fi;
done < $horse_trans/working_list_NoPBMCs_NoCereb.txt

## create list of assemblies for tissues of multiple libraries
rm -f $tissue_Digimerge/tissue_assemblies.txt
for tissue in $tissue_Cuffmerge/*; do
  echo ${tissue#$tissue_Digimerge/}/withoutGTF >> $tissue_Digimerge/digi_tissue_assemblies.txt;
  echo ${tissue}/withoutGTF >> $horse_trans/rawdigi_TopCuff_assemblies.txt; done

## initiate a given track hub
hub_name=$"HorseTrans3"
shortlabel=$"rawdigi_TopCuff"
longlabel=$"diginorm of raw data followed by reference guided Tophat/Cufflinks"
email=$"drtamermansour@gmail.com"
cd $track_hub
bash $script_path/create_trackHub.sh "$UCSCgenome" "$hub_name" "$shortlabel" "$longlabel" "$email"

### Run Cufflinks: output transcripts.gtf in the same tophat_sample folder
#while read work_dir; do
#  echo $work_dir
#  cd $work_dir/digi_tophat_output
#  sample_list=$work_dir/digi_tophat_output/sample_list.txt
#  bash ${script_path}/run_cufflinks.sh "$sample_list" "$refGTF_file" "$script_path/cufflinks2.sh";
#done < $horse_trans/working_list_NoPBMCs.txt

## Check for successful Cufflinks runs and trouble shooting the failed Cufflinks jobs (requires cufflinks.e)
#while read work_dir; do
#  cd $work_dir/digi_tophat_output
#  sample_list=$work_dir/digi_tophat_output/failedSamples.txt           ## define the path of empty file
#  bash $script_path/check_cufflinks.sh "$sample_list"
#  x=$(cat $sample_list | wc -l)
#  if [ $x -ne 0 ]; then
#    echo "Failed Cufflinks jobs in: "$work_dir
#    cat $sample_list
#     bash ${script_path}/run_cufflinks.sh "$sample_list" "$refGTF_file" "$script_path/cufflinks.sh";
#fi; done < $horse_trans/working_list_NoPBMCs.txt
###################
## Assess computational utilization of cufflinks
#cufflinks_utlization=$prepData/digi_cufflinks_utlization.txt
#> $cufflinks_utlization
#while read work_dir; do
#  cd $work_dir/digi_tophat_output
#  echo $(pwd) >> $cufflinks_utlization
#  for dir in tophat_*; do
#    f=$dir/cufflinks.e*
#    if [ -f $f ]; then
#      f2=$(echo $f | sed 's/cufflinks.e/cufflinks.o/')
#      echo $dir >> $cufflinks_utlization
#      grep "resources_used.vmem =" $f2 >> $cufflinks_utlization
#      grep "resources_used.walltime =" $f2 >> $cufflinks_utlization
#      grep "Resource_List.nodes =" $f2 >> $cufflinks_utlization
#    else echo "no cufflinks report in $dir" >> $cufflinks_utlization
#fi; done; done < $horse_trans/working_list_NoPBMCs.txt

###########################################################################################

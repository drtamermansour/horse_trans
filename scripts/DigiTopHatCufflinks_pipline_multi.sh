#!/bin/sh
myRoot=$"/mnt/ls12/Tamer"
source $myRoot/config.txt
###########################################################################################
## pipeline_diginormAllsamples_mergeSamples_Tophat2.nonGuided_Cufflinks
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
#cutoff=50
cutoff=10
normReads="normalizied_RNA_reads_k${kmer}_C${cutoff}"
while read work_dir; do
  echo $work_dir
  mkdir -p $work_dir/$normReads
  cd $work_dir/$normReads
  lib=$(basename $work_dir | cut -d"_" -f 1)
  bash ${script_path}/run_diginorm.sh "$lib" "$work_dir/trimmed_RNA_reads" "$kmer" "$cutoff" "$script_path"
done < $horse_trans/working_list_NoPBMCs.txt

## Check for successful diginorm and trouble shooting the failed jobs (requires T_Trim.e)
sample_list=$prepData/failed_diginorm_k${kmer}_C${cutoff}.txt            ## define the path of empty file
> $sample_list
while read work_dir; do
  cd $work_dir/$normReads
  bash $script_path/check_diginorm.sh "$sample_list"
done < $horse_trans/working_list_NoPBMCs.txt
x=$(cat $sample_list | wc -l)
if [ $x -ne 0 ]; then
  echo "Failed jobs in: "
  cat $sample_list
  while read work_dir; do
    cd $work_dir/$normReads
    lib=$(basename $work_dir | cut -d"_" -f 1)
    bash ${script_path}/run_diginorm.sh "$lib" "$work_dir/trimmed_RNA_reads" "$kmer" "$cutoff" "$script_path"
  done < $sample_list
fi

## split the interleaved reads
while read work_dir; do
  echo $work_dir
  cd $work_dir/$normReads
  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
  if [ "$lib" = $"PE" ]; then
    sample_list=$work_dir/trimmed_RNA_reads/sample_list.txt
    bash ${script_path}/run_split_reads.sh "$sample_list" $script_path/split_reads.sh
fi; done < $horse_trans/working_list_NoPBMCs.txt

## merge singletons and change the file names to fit the tophat script
#while read work_dir; do
#  echo $work_dir
#  cd $work_dir/$normReads
#  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
#  if [ "$lib" = $"PE" ]; then
#    singletones=1
#    for f in *_R1_001.pe.fq; do
#      base=${f%_R1_001.pe.fq}
#      if [ $singletones -eq 1 ]; then cat $f allsingletons.fq.keep > "$base"_R1_001.pe.se.fq; singletones=0;
#      else mv $f "$base"_R1_001.pe.se.fq; fi; done
#  elif [ "$lib" = $"SE" ]; then
#    mv allsingletons.fq.keep allsingletons_SR_002.se.fq
#fi; done < $horse_trans/working_list_NoPBMCs.txt

## merge the files in tophat compatible format
while read work_dir; do
  echo $work_dir
  cd $work_dir/$normReads
  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
  if [ "$lib" = $"PE" ]; then
    cat *_R1_001.pe.fq allsingletons.fq.keep > allsamples_R1_002.pe.se.fq
    cat *_R2_001.pe.fq > allsamples_R2_002.pe.fq
  elif [ "$lib" = $"SE" ]; then
    mv allsingletons.fq.keep allsingletons_SR_002.se.fq
fi; done < $horse_trans/working_list_NoPBMCs.txt

## merge the files in tophat compatible format
#while read work_dir; do
#  echo $work_dir
#  mkdir -p $work_dir/${normReads}.merged
#  cd $work_dir/${normReads}.merged
#  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
#  if [ "$lib" = $"PE" ]; then
#    cat *_R1_001.pe.fq >> allsamples_R1_002.pe.fq
#    cat *_R2_001.pe.fq >> allsamples_R2_002.pe.fq
#    cat allsingletons.fq.keep >> allsingletons_R_002.se.fq
#  elif [ "$lib" = $"SE" ]; then
#    cat allsingletons.fq.keep >> allsingletons_R_002.se.fq
#fi; done < $horse_trans/working_list_NoPBMCs.txt


## define the list samples.
## This is where you can edit the output list file(s) to restrict the processing for certain target(s)
while read work_dir; do if [ -d $work_dir/$normReads ]; then
  rm -f $work_dir/$normReads/sample_list.txt
  for f in $work_dir/$normReads/{*_R1_*.pe.se.fq,*_SR_*.se.fq}; do if [ -f $f ]; then
    echo $f >> $work_dir/$normReads/sample_list.txt; fi; done;
fi; done < $horse_trans/working_list_NoPBMCs.txt

## run Tophat on each library as a one sample

while read work_dir; do
  echo $work_dir
  mkdir -p $work_dir/digiMulti.k${kmer}.C${cutoff}_tophat_output
  cd $work_dir/digiMulti.k${kmer}.C${cutoff}_tophat_output
  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
  strand=$(basename $work_dir | cut -d"_" -f 3 | sed 's/\./-/')   ## fr-unstranded, fr-firststrand or fr-secondstrand
  sample_list=$work_dir/$normReads/sample_list.txt
  bash ${script_path}/run_tophat.sh "$sample_list" "$lib" "$strand" "$Bowtie2_genome_index_base" "$transcriptome_index" "$script_path"
done < $horse_trans/working_list_NoPBMCs.txt

## Check for successful tophat runs and trouble shooting the failed tophat jobs (require tophat-[SP]E.e & .o)
while read work_dir; do
  cd $work_dir/digiMulti.k${kmer}.C${cutoff}_tophat_output
  sample_list=$work_dir/trimmed_RNA_reads/sample_list.txt
  failedSample_list=$work_dir/trimmed_RNA_reads/tophat_failedSamples.txt
  > $failedSample_list                                                ## erase previouslly failed samples if any
  bash $script_path/check_tophat.sh "$failedSample_list"        ## require tophat-[SP]E.e & .o
  bash $script_path/check2_tophat.sh "$sample_list" "$failedSample_list"   ## check output log files
  x=$(cat $failedSample_list | wc -l)
  if [ $x -ne 0 ]; then
    lib=$(basename $work_dir | cut -d"_" -f 1)
    strand=$(basename $work_dir | cut -d"_" -f 3 | sed 's/\./-/')
    echo "Failed tophat jobs in: "$work_dir
#    bash ${script_path}/run_tophat.sh "$failedSample_list" "$lib" "$strand" "$Bowtie2_genome_index_base" "$transcriptome_index" "$script_path"
fi; done < $horse_trans/working_list_NoPBMCs.txt

##################
## create summary for tophat run
headers=$(Rscript -e 'cat("Tissue", "Library", "min_mapping", "max_mapping", "min_concordance", "max_concordance", sep="\t");')
echo "$headers" > $horse_trans/digiMulti.k${kmer}.C${cutoff}_tophat_summary.txt
while read work_dir; do
  > $work_dir/digiMulti.k${kmer}.C${cutoff}_tophat_output/allsample_summary.txt
  for f in $work_dir/digiMulti.k${kmer}.C${cutoff}_tophat_output/tophat_*; do
    echo ${f} >> $work_dir/digiMulti.k${kmer}.C${cutoff}_tophat_output/allsample_summary.txt
    cd ${f}
    grep "overall read mapping rate" align_summary.txt >> ../allsample_summary.txt
    grep "concordant pair alignment rate" align_summary.txt >> ../allsample_summary.txt
  done
  mapping=$(grep "overall read mapping rate" $work_dir/digiMulti.k${kmer}.C${cutoff}_tophat_output/allsample_summary.txt | awk '{ print $1 }' | sort -n | sed -e 1b -e '$!d' | tr "\n" "\t")
  conc=$(grep "concordant pair alignment rate" $work_dir/digiMulti.k${kmer}.C${cutoff}_tophat_output/allsample_summary.txt | awk '{ print $1 }' | sort -n | sed -e 1b -e '$!d' | tr "\n" "\t")
  lib=$(basename $work_dir)
  tissue=$(dirname $work_dir | xargs basename)
  echo "$tissue"$'\t'"$lib"$'\t'"$mapping""$conc" >> $horse_trans/digiMulti.k${kmer}.C${cutoff}_tophat_summary.txt
done < $horse_trans/working_list_NoPBMCs.txt

##################
## define the list samples.
## This is where you can edit the output list file(s) to restrict the processing for certain target(s)
rm -f $horse_trans/digiMulti.k${kmer}.C${cutoff}_sample_list.txt
while read work_dir; do if [ -d $work_dir/digiMulti.k${kmer}.C${cutoff}_tophat_output ]; then
  for f in $work_dir/digiMulti.k${kmer}.C${cutoff}_tophat_output/tophat_*; do if [ -d $f ]; then
    echo $f >> $horse_trans/digiMulti.k${kmer}.C${cutoff}_sample_list.txt; fi; done;
fi; done < $horse_trans/working_list_NoPBMCs.txt

## Merge BAM files
tissue_Digimerge=$tissue_merge/digimerge
mkdir -p $tissue_Digimerge/digiMulti.k${kmer}.C${cutoff}
cd $tissue_Digimerge/digiMulti.k${kmer}.C${cutoff}
samples=()
while read sample; do samples+=($sample/accepted_hits.bam); done < $horse_trans/digiMulti.k${kmer}.C${cutoff}_sample_list.txt
len=${#samples[@]}
if [ $len -gt 1 ]; then module load SAMTools/0.1.19; samtools merge merged.bam ${samples[*]};
elif [ $len -eq 1 ]; then cp ${samples[0]} merged.bam;
else echo "can find bam files"; fi




### Run Cufflinks: output transcripts.gtf in the same tophat_sample folder
cufflinks_run="nonGuided_Cufflinks"
cd $tissue_Digimerge/digiMulti.k${kmer}.C${cutoff}
sample=$tissue_Digimerge/digiMulti.k${kmer}.C${cutoff}/merged.bam
label="digiMulti.k${kmer}.C${cutoff}"
bash ${script_path}/run_cufflinks_noRef_single.sh "$sample" "$label" "$script_path/cufflinks_noRef2.sh";




## Check for successful Cufflinks runs and trouble shooting the failed Cufflinks jobs (requires cufflinks.e)
cd $tissue_Digimerge/digiMulti.k${kmer}.C${cutoff}
failedSample_list=$horse_trans/digiMulti.k${kmer}.C${cutoff}_failedSamples.txt           ## define the path of empty file
bash $script_path/check_cufflinks2.sh "$failedSample_list"
x=$(cat $failedSample_list | wc -l)
if [ $x -ne 0 ]; then
  echo "Failed Cufflinks jobs in: "$horse_trans/digiMulti.k${kmer}.C${cutoff}
  cat $failedSample_list
  bash ${script_path}/run_cufflinks_noRef_single.sh "$sample" "$label" "$script_path/cufflinks_noRef2.sh";
fi
##################
## relocate the cufflinks analysis results
cd $tissue_Digimerge/digiMulti.k${kmer}.C${cutoff}
mkdir $cufflinks_run && \
mv $tissue_Digimerge/digiMulti.k${kmer}.C${cutoff}/{transcripts.gtf,skipped.gtf,*.fpkm_tracking,cufflinks.[oe]*} $cufflinks_run/.
##################
## Add to the list of assemblies for tissues of multiple libraries
for tissue in $tissue_Digimerge/digiMulti.k${kmer}.C${cutoff}/$cufflinks_run; do if [ -d $tissue ]; then
  echo ${tissue#$tissue_Digimerge/} >> $tissue_Cuffmerge/${cufflinks_run}_${cuffmerge_run}_tissue_assemblies.txt;
fi; done
####################
## convert the gtf files into BigBed files & copy the BigBed files to the track hub directory
assembly="digiMulti.k20.C10/nonGuided_Cufflinks"
cd $tissue_Digimerge/$assembly
targetAss=$"transcripts.gtf"
bash $script_path/gtfToBigBed.sh "$targetAss" "$genome_dir/$UCSCgenome.chrom.sizes" "$script_path"
if [ -f $"transcripts.BigBed" ];then
  identifier=$(echo $assembly | sed 's/\//_/g')
  cp transcripts.BigBed $track_hub/$UCSCgenome/BigBed/${identifier}.BigBed
  echo $tissue_Digimerge/$assembly >> $horse_trans/Tophat_${cufflinks_run}_${cuffmerge_run}_assemblies.txt;
fi

## Add to the HorseTrans_TopNonGuidCuff track hub
shortlabel=$"TopNonGuidCuff"
current_libs=$track_hub/current_libs_$shortlabel
current_tissues=$track_hub/current_tiss_$shortlabel
trackDb=$track_hub/$UCSCgenome/trackDb_$shortlabel.txt
lib_assemblies=$prepData/${cufflinks_run}_${cuffmerge_run}_merged_assemblies.txt
tiss_assemblies=$tissue_Cuffmerge/${cufflinks_run}_${cuffmerge_run}_tissue_assemblies.txt
bash $script_path/edit_trackDb.sh $current_libs $current_tissues $trackDb $lib_assemblies $tiss_assemblies

## create the HTML file page for every track
###########################################################################################

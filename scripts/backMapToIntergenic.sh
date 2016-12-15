cd horse_trans
horse_trans=$(pwd)
source $horse_trans/user_config.txt
source $horse_trans/config.txt
tissue_Cuffmerge=$tissue_merge/cuffmerge
dist_dir="all_tissues_frac$isoformfrac"
cuffmerge_output=$dist_dir/$cufflinks_run/$cuffmerge_run
filtered2_hiExpGTF=$tissue_Cuffmerge/$cuffmerge_output/filtered/highExp/merged.gtf

## BAck mapping to the intergenic transcripts
## make a version of the refined transcriptome with articial version of mitochondrial sequences 
mkdir $horse_trans/backMap_intergenic
cd $horse_trans/backMap_intergenic
grep -F -w -f $tissue_Cuffmerge/$cuffmerge_output/filtered/supported/unsup $filtered2_hiExpGTF > merged.gtf ## 38507
assembly="$horse_trans/backMap_intergenic/merged.gtf"
bash $script_path/run_genome_to_cdna_fasta.sh "merged.gtf" "$genome" "transcripts.fasta" "$script_path/genome_to_cdna_fasta.sh"
bash $script_path/run_salmonIndex.sh "horse_index" "transcripts.fasta" ${script_path}/salmonIndex.sh
while read work_dir; do
  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
  strand=$(basename $work_dir | cut -d"_" -f 3 | sed 's/\./-/')   ## fr-unstranded, fr-firststrand or fr-secondstrand
  identifier=$(echo $work_dir | rev | cut -d"/" -f 1,2 | rev | sed 's/\//_/')
  seq_dir=$work_dir/trimmed_RNA_reads
  bash ${script_path}/run_salmon.sh "$lib" "$strand" "horse_index" "$identifier" "$seq_dir" "$script_path"
done < $horse_trans/working_list.txt
find ./*.quant -name *.sf -exec grep -H "mapping rate" {} \; | sort > salmonQuant_summary.txt
echo "Total no of input reads" > salmonQuant_summary_detailed.txt
find ./*.quant -name salmon_quant.log -exec grep -H "total fragments" {} \; | sort >> salmonQuant_summary_detailed.txt
echo "No of mapped reads" >> salmonQuant_summary_detailed.txt
find ./*.quant -name salmon_quant.log -exec grep -H "total reads" {} \; | sort >> salmonQuant_summary_detailed.txt
echo "Mapping rate" >> salmonQuant_summary_detailed.txt
find ./*.quant -name salmon_quant.log -exec grep -H "Mapping rate" {} \; | sort >> salmonQuant_summary_detailed.txt
#cp salmonQuant_summary*.txt $horse_trans/downloads/backmapping_stats/.
python $script_path/gather-counts2.py -i "$(pwd)"
echo "transcript"$'\t'"length" > transcripts.lengthes
sf=$(find ./*.quant -name \*.sf | head -n1)
cat $sf | grep -v "^#" | awk -F "\t" -v OFS='\t' '{print $1,$2}' >> transcripts.lengthes
grep "^>" transcripts.fasta | awk -F'[> ]' '{print $3,$2}' > gene_transcript_map

## library specific expression and assembly
#echo "## no of genes and transcripts" > $horse_trans/libAsmStats.txt
while read work_dir; do
  tissue=$(basename $(dirname $work_dir))
  lib=$(basename $work_dir)
  target=$tissue"_"$lib
  echo $target
  bash $script_path/run_calcTPM.sh "$(pwd)" "$target" "transcripts.lengthes" "gene_transcript_map" ${script_path}/calcTPM2.R
#  dir=$work_dir/tophat_output/$cufflinks_run/$cuffmerge_run/filtered
#  mkdir -p $dir
#  cat $target.dataSummary_comp | tail -n+2 | awk '{if($10 >= 5)print $3}' > $target.keepit.id
#  grep -F -w -f $target.keepit.id $assembly | awk '$1 != "chrM"' > $dir/merged.gtf
#  ## copy the annotation to the download folder
#  cp $dir/merged.gtf $horse_trans/downloads/tissueSpecific_assemblies/$target.gtf
#  ## statistics (no of genes and transcripts)
#  cat $dir/merged.gtf | awk -F '[\t"]' '{ print $10 }' |  sort | uniq | wc -l
#  cat $dir/merged.gtf | awk -F '[\t"]' '{ print $12 }' |  sort | uniq | wc -l
#  ## self cuffcompare
#  bash ${script_path}/run_cuffcompareSelf.sh "$dir/merged.gtf" "selfCuff_$target" "selfCuff_$target.log" "$script_path/cuffcompareSelf.sh"
done < $horse_trans/working_list.txt #>> $horse_trans/libAsmStats.txt

## Tissue specific expression and assembly
#echo "## no of genes and transcripts" > $horse_trans/tisAsmStats.txt
while read work_dir; do
  target=$(basename $work_dir)
  echo $target
  bash $script_path/run_calcTPM.sh "$(pwd)" "$target" "transcripts.lengthes" "gene_transcript_map" ${script_path}/calcTPM2.R
#  dir=$tissue_Cuffmerge/$target/$cufflinks_run/$cuffmerge_run/filtered
#  mkdir -p $dir
#  cat $target.dataSummary_comp | tail -n+2 | awk '{if($10 >= 5)print $3}' > $target.keepit.id
#  grep -F -w -f $target.keepit.id $assembly | awk '$1 != "chrM"' > $dir/merged.gtf
#  ## copy the annotation to the download folder
#  cp $dir/merged.gtf $horse_trans/downloads/tissueSpecific_assemblies/$target.gtf
#  ## statistics
#  cat $dir/merged.gtf | awk -F '[\t"]' '{ print $10 }' |  sort | uniq | wc -l
#  cat $dir/merged.gtf | awk -F '[\t"]' '{ print $12 }' |  sort | uniq | wc -l
#  ## self cuffcompare
#  bash ${script_path}/run_cuffcompareSelf.sh "$dir/merged.gtf" "self" "cuffcompare_self.log" "$script_path/cuffcompareSelf.sh"
done < $horse_trans/multi_lib_tissues.txt #>> $horse_trans/tisAsmStats.txt
#cp $horse_trans/libAsmStats.txt $horse_trans/tisAsmStats.txt $horse_trans/downloads/backmapping_stats/.

## calculate tissue specific expression
targets=()
i=1
rm temp.* isotemp.*
while read target;do
  f="$target".dataSummary_comp
  if [ ! -f $f ];then f=$(ls "$target"_*.dataSummary_comp);fi
  cat $f | tail -n+2 | awk '{print $2,$7}' | uniq > $f.gene
  cat $f | tail -n+2 | awk '{print $3,$6}' > $f.isoform
  targets+=($target)
  if [ $i -eq 1 ];then cat $f.gene > temp.$i;else join -t" " --nocheck-order temp.$((i-1)) $f.gene > temp.$i;fi
  if [ $i -eq 1 ];then cat $f.isoform > isotemp.$i;else join -t" " --nocheck-order isotemp.$((i-1)) $f.isoform > isotemp.$i;fi
  ((i+=1))
done < <(ls *_*.dataSummary_comp | awk -F '_' '{print $1}' | sort | uniq)
echo "geneName" "${targets[@]}" > allTissues_geneTPM
cat temp.$((i-1)) >> allTissues_geneTPM
echo "isoformName" "${targets[@]}" > allTissues_isoformTPM
cat isotemp.$((i-1)) >> allTissues_isoformTPM
rm temp.* isotemp.*


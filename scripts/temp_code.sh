## Get all accepted hits
samtools view accepted_hits.bam > allaccepted.sam
cat allaccepted.sam | wc -l	#4917873
## Get them as , 1st, 2nd and unpaired reads 
samtools view -f 64 accepted_hits.bam > firstpair.sam
cat firstpair.sam | wc -l	#2562200
samtools view -f 128 accepted_hits.bam > secondpair.sam
cat secondpair.sam | wc -l	#2267793
samtools view -F 1 accepted_hits.bam > unpaired.sam
cat unpaired.sam | wc -l	#87880


cat firstpair.sam | cut -f1 | sort | uniq | wc -l	#2500346 / 22314202
cat secondpair.sam | cut -f1 | sort | uniq | wc -l	#2200986 / 22314202
cat unpaired.sam | cut -f1 | sort | uniq | wc -l	#85389 / 978030

samtools view -F 256 -f 64 accepted_hits.bam > firstpair_1ry.sam
cat firstpair_1ry.sam | wc -l	#2500346


###########################################################################################

## Try the extra senstive Bowtie mapping

script_path=$myRoot/horse_trans/scripts2


## define the list samples.
## This is where you can edit the output list file(s) to restrict the processing for certain target(s)
while read work_dir; do if [ -d $work_dir/trimmed_RNA_reads ]; then
  rm -f $work_dir/trimmed_RNA_reads/sample_list.txt
  for f in $work_dir/trimmed_RNA_reads/{*_R1_*.pe.se.fq,*_SR_*.se.fq}; do if [ -f $f ]; then
    echo $f >> $work_dir/trimmed_RNA_reads/sample_list.txt; fi; done;
  fi; done < $horse_trans/working_list3.txt


## run Tophat on each sample
while read work_dir; do
  echo $work_dir
  mkdir -p $work_dir/tophat_output_sens
  cd $work_dir/tophat_output_sens
  lib=$(basename $work_dir | cut -d"_" -f 1)                      ## PE or SE
  strand=$(basename $work_dir | cut -d"_" -f 3 | sed 's/\./-/')   ## fr-unstranded, fr-firststrand or fr-secondstrand
  sample_list=$work_dir/trimmed_RNA_reads/sample_list.txt
  bash ${script_path}/run_tophat.sh "$sample_list" "$lib" "$strand" "$Bowtie2_genome_index_base" "$transcriptome_index" "$script_path"
done < $horse_trans/working_list3.txt

#################################################################################
##### Suggested piplines
### With Tophat
## Compare Ref-guided versus reference free
## Try very senstive Bowtie search
## Try --coverage-search
## For diginorm, try to start with mapped reads only (to minimize error accumulation ==> decrease RAM required to achive minimal false postive rate and allow either final mapping). You can diginorm the non-mapped reads separetly



###########################################################################################

## Split GTF and BAM files by chromosomes
cat ${refGTF_file} | awk '{ print $1}' | sort | uniq > $genome_dir/refChrom_list.txt
mkdir $genome_dir/GTFbyChrom
while read refChrom; do
  cat ${refGTF_file} | awk -v ref="$refChrom" '($1 == ref)' >> $genome_dir/GTFbyChrom/"$refChrom"_refGene.gtf
done < $genome_dir/refChrom_list.txt

module load BAMTools/2.2.3
bamtools split -in accepted_hits.bam -reference

module load cufflinks/2.2.1
cufflinks --GTF-guide $genome_dir/GTFbyChrom/chr10_refGene.gtf --num-threads 4 --verbose accepted_hits.REF_chr10.bam















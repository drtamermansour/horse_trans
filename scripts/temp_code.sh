myRoot=$"/mnt/ls12/Tamer"
data_path=$myRoot/horse_trans/data2
work_dir=$data_path/PE_UNS_Bellone_Retina
genome_dir=$myRoot/horse_trans/refGenome
Bowtie2_genome_index_base=$genome_dir/Bowtie2Index/genome
transcriptome_index=$genome_dir/trans_index/equ
script_path=$myRoot/horse_trans/scripts/temp
lib_type=$"fr-unstranded"

f1=${work_dir}/trimmed_RNA_reads/*_R1_*.pe.se.fq
f2=${work_dir}/trimmed_RNA_reads/*_R2_*.pe.fq
output_dir=$work_dir/tophat_output/repeat2
mkdir ${output_dir}
qsub -v output_dir="${output_dir}",\
Bowtie2_genome_index_base="${Bowtie2_genome_index_base}",\
transcriptome_index="${transcriptome_index}",\
lib_type="${lib_type}",\
f1="$f1",f2="$f2" $script_path/tophat.sh


f1=${work_dir}/trimmed_RNA_reads2/*_R1_*.pe.fq
f2=${work_dir}/trimmed_RNA_reads2/*_R2_*.pe.fq
output_dir=$work_dir/tophat_output/noSingletons
mkdir ${output_dir}
qsub -v output_dir="${output_dir}",\
Bowtie2_genome_index_base="${Bowtie2_genome_index_base}",\
transcriptome_index="${transcriptome_index}",\
lib_type="${lib_type}",\
f1="$f1",f2="$f2" $script_path/tophat2.sh

f1=${work_dir}/trimmed_RNA_reads4/*_R1_*.pe.se.fq
f2=${work_dir}/trimmed_RNA_reads4/*_R2_*.pe.fq
output_dir=$work_dir/tophat_output/newVersion
mkdir ${output_dir}
qsub -v output_dir="${output_dir}",\
Bowtie2_genome_index_base="${Bowtie2_genome_index_base}",\
transcriptome_index="${transcriptome_index}",\
lib_type="${lib_type}",\
f1="$f1",f2="$f2" $script_path/tophat4.sh

f1=${work_dir}/trimmed_RNA_reads5/*_R1_*.pe.se.fq
f2=${work_dir}/trimmed_RNA_reads5/*_R2_*.pe.fq
output_dir=$work_dir/tophat_output/simplified
mkdir ${output_dir}
qsub -v output_dir="${output_dir}",\
Bowtie2_genome_index_base="${Bowtie2_genome_index_base}",\
transcriptome_index="${transcriptome_index}",\
lib_type="${lib_type}",\
f1="$f1",f2="$f2" $script_path/tophat5.sh

#####
f1=${work_dir}/trimmed_RNA_reads3/*_R1_*.pe.se.fq
f2=${work_dir}/trimmed_RNA_reads3/*_R2_*.pe.fq
output_dir=$work_dir/tophat_output/noTranscriptome
mkdir ${output_dir}
qsub -v output_dir="${output_dir}",\
Bowtie2_genome_index_base="${Bowtie2_genome_index_base}",\
transcriptome_index="${transcriptome_index}",\
lib_type="${lib_type}",\
f1="$f1",f2="$f2" $script_path/tophat3.sh

f1=${work_dir}/trimmed_RNA_reads6/*_R1_*.pe.se.fq
f2=${work_dir}/trimmed_RNA_reads6/*_R2_*.pe.fq
output_dir=$work_dir/tophat_output/noCovsearch
mkdir ${output_dir}
qsub -v output_dir="${output_dir}",\
Bowtie2_genome_index_base="${Bowtie2_genome_index_base}",\
transcriptome_index="${transcriptome_index}",\
lib_type="${lib_type}",\
f1="$f1",f2="$f2" $script_path/tophat6.sh

f1=${work_dir}/trimmed_RNA_reads7/*_R1_*.pe.se.fq
f2=${work_dir}/trimmed_RNA_reads7/*_R2_*.pe.fq
output_dir=$work_dir/tophat_output/noMicrosearch
mkdir ${output_dir}
qsub -v output_dir="${output_dir}",\
Bowtie2_genome_index_base="${Bowtie2_genome_index_base}",\
transcriptome_index="${transcriptome_index}",\
lib_type="${lib_type}",\
f1="$f1",f2="$f2" $script_path/tophat7.sh

f1=${work_dir}/trimmed_RNA_reads8/*_R1_*.pe.se.fq
f2=${work_dir}/trimmed_RNA_reads8/*_R2_*.pe.fq
output_dir=$work_dir/tophat_output/noPrefilter
mkdir ${output_dir}
qsub -v output_dir="${output_dir}",\
Bowtie2_genome_index_base="${Bowtie2_genome_index_base}",\
transcriptome_index="${transcriptome_index}",\
lib_type="${lib_type}",\
f1="$f1",f2="$f2" $script_path/tophat8.sh
#############################################################

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


###################################################################

cd ${work_dir}/b_diginorm$suffix
cat *.keep > allsamples.keep.fq
qsub -v input_file=$"allsamples.keep.fq" ${script_path}/jellyfish.sh

## filter abundance
norm_files=()
for f in *.keep; do norm_files+=($f); done;
qsub -v input=$"allsamples$suffix.kt",files="${norm_files[*]}" ${script_path}/filter_abund.sh
#filter-abund.py -V $input $files
for f in filter_abund.e*; do grep -B 2 "^output in" $f >> filter_abund.summary; done

## break out the orphaned and still-paired reads & rename files (this step ends with .s_pe.fq & .s_se.fq for each sample)
for i in J*.s_pe.*.abundfilt; do extract-paired-reads.py $i; done
for i in P_ast*.s_pe.*.abundfilt; do qsub -v input=$i $script_path/extract-paired-reads.sh; done
##  combine the orphaned reads into a single file & rename pe files
for i in *.s_se.fq.keep.abundfilt; do
pe_orphans=$(basename $i .s_se.fq.keep.abundfilt).s_pe.fq.keep.abundfilt.se
cat $i $pe_orphans > $(basename $i .keep.abundfilt)
done
rm *.abundfilt.se
for i in *.abundfilt.pe; do mv $i $(basename $i .keep.abundfilt.pe); done
## split interleaved files & ## merge the single reads to the end of left reads (this step reform the data into .s_pe1_se.fq & .s_pe2.fq for each sample )
for i in *.s_pe.fq; do split-paired-reads.py $i; done
for f in *.s_pe.fq; do
echo $f  >> check_file_ends;
base=$(basename $f);
tail -n 8 $f | head -n 4 >> check_file_ends;
tail -n 4 $base.1 >> check_file_ends;
tail -n 4 $f >> check_file_ends;
tail -n 4 $base.2 >> check_file_ends;
done
for f in *.s_pe.fq.1; do echo $f; base=$(basename $f .s_pe.fq.1); cat $f $base.s_se.fq > $base.s_pe1_se.fq; done
#rm *.s_pe.fq.1
for f in *.s_pe.fq.2; do mv $f $(basename $f .s_pe.fq.2).s_pe2.fq; done

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

















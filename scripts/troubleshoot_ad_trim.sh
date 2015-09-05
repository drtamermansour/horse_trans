work_dir="$1"
script_path="$2"

## Troubleshooting
cd $work_dir/trimmed_RNA_reads/job_reports
while read line; do
R1=$(grep "Started with arguments" $line | cut -d" " -f 8)
basename $R1 | cut -d"_" -f1 >> failedSamples
done < failedJobs

cd ../
mkdir -p failedSamples
while read sample; do
mv "$sample"_* failedSamples/.
done < job_reports/failedSamples

while read sample; do
f=$"$work_dir/fastq_data/$sample"*_R1_*.fastq.gz
input_one="$f"
echo $input_one
input_two=$(echo "$f" | sed s/_R1_/_R2_/)
#echo $input_two

temp=$(basename $f)
temp2=${temp%.fastq.gz}

output_pe1=$work_dir/trimmed_RNA_reads/$temp2".pe.fq"
#echo $output_pe1
output_pe2=$(echo "$output_pe1" | sed s/_R1_/_R2_/)
#echo $output_pe2
output_se1=$work_dir/trimmed_RNA_reads/$temp2".se.fq"
#echo $output_se1
output_se2=$(echo "$output_se1" | sed s/_R1_/_R2_/)
#echo $output_se2
qsub -v R1_INPUT="$input_one",R2_INPUT="$input_two",output_pe1="$output_pe1",output_pe2="$output_pe2",output_se1="$output_se1",output_se2="$output_se2" ${script_path}/adapter_trimmer.sh
done < $work_dir/trimmed_RNA_reads/job_reports/failedSamples

## change /2 to /1 in s2_se then combine single reads
while read sample; do
f=$"$work_dir/trimmed_RNA_reads/$sample"*_R1_*.se.fq
echo $f
f2=$(basename $f | sed 's/_R1_/_R2_/');
newf=$(basename $f | sed 's/_R1_/_R_/');
sed 's/\/2$/\/1/g' $f2 > $f2.temp;
cat $f $f2.temp > $newf;
rm *.se.fq.temp
done < $work_dir/trimmed_RNA_reads/job_reports/failedSamples

## merge the single reads to the end of s1_pe file to make s1_pe_se
while read sample; do
f=$"$work_dir/trimmed_RNA_reads/$sample"*_R1_*.pe.fq
fr=$(basename $f | sed 's/_R1_/_R_/');
fr2=$(echo $fr | sed 's/.pe.fq/.se.fq/');
newf=$(basename $f | sed 's/.pe.fq/.pe.se.fq/');
cat $f $fr2 > $newf;
done < $work_dir/trimmed_RNA_reads/job_reports/failedSamples

## remove failed tophat jobs
mkdir $work_dir/tophat_output/failedTopHatSamples
while read sample; do
output_dir=$work_dir/tophat_output"/tophat_"$sample
echo $output_dir
mv $output_dir $work_dir/tophat_output/failedTopHatSamples/.
done < $work_dir/trimmed_RNA_reads/job_reports/failedSamples

## run tophat for missed samples
myRoot=$"/mnt/ls12/Tamer"
genome_dir=$myRoot/horse_trans/refGenome
Bowtie2_genome_index_base=$genome_dir/Bowtie2Index/genome
transcriptome_index=$genome_dir/trans_index/equ
lib_type=$"fr-unstranded"
while read sample; do
f1=$"$work_dir/trimmed_RNA_reads/$sample"*_R1_*.pe.se.fq
f2=$"$work_dir/trimmed_RNA_reads/$sample"*_R2_*.pe.fq
base=$(basename $f1)
base2=${base%_R1_*.pe.se.fq}
output_dir=$work_dir/tophat_output"/tophat_"$base2
echo $output_dir
mkdir ${output_dir}
qsub -v output_dir="${output_dir}",\
Bowtie2_genome_index_base="${Bowtie2_genome_index_base}",\
transcriptome_index="${transcriptome_index}",\
lib_type="${lib_type}",\
f1=$f1,f2=$f2 $script_path/tophat.sh
done < $work_dir/trimmed_RNA_reads/job_reports/failedSamples








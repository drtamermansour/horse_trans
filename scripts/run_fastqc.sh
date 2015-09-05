RNA_reads="$1"
# script="$2"

## QC of fastq data
module load FastQC
cd $RNA_reads
for f in $RNA_reads/*.fastq.gz; do
    echo $f >> report
    fastqc -f fastq -noextract $f
    #qsub -v INPUT_FILE="$f" $script
    html_temp1=$(basename $f)
    html_temp2=${html_temp1%.fastq.gz}
    html_file=$html_temp2"_fastqc.html"
    tr '>' '\n' < $html_file > temp
    tr '<' '\n' < temp > temp2
    encode=$(grep -A 4 "Encoding" temp2 | head -5| tail -1)
    echo $encode >> report
    total_seq=$(grep -A 4 "Total Sequences" temp2 | head -5 | tail -1)
    echo $total_seq >> report
    seq_len=$(grep -A 4 "Sequence length" temp2 | head -5 | tail -1)
    echo $seq_len >> report
    rm temp*
done

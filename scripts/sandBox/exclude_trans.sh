trans="$1"
label="$2"
gene_list="$3"
module load QIIME/1.8.0
filter_fasta.py --input_fasta_fp $trans --output_fasta_fp $label.fasta --seq_id_fp $gene_list  --negate


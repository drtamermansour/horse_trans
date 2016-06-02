R1_INPUT="$1"
R2_INPUT="$2"
output_pe1="$3"
output_pe2="$4"
output_se1="$5"
output_se2="$6"

script_path=$(dirname "${BASH_SOURCE[0]}")
TrimmomaticPE -threads 4 -phred33 ${R1_INPUT} ${R2_INPUT} ${output_pe1} ${output_se1} ${output_pe2} ${output_se2} ILLUMINACLIP:$script_path/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:2 MINLEN:20 &> ${R1_INPUT}.log

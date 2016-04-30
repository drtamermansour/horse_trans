work_path="$1"
identifier="$2" 
transLen="$3" 
gene_transcript_map="$4"
script="$5"

module load R/3.0.1

Rscript "$script" "$work_path" "$identifier" "$transLen" "$gene_transcript_map" > _tmp


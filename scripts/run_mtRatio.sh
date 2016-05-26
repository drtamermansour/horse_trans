geneTPM="$1"
output="$2"
script="$3"

module load R/3.2.0

Rscript $script $geneTPM $output


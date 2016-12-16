geneTPM="$1"
output="$2"
script="$3"

export R_LIBS_USER=~/R/v3.2.0/library
module load R/3.2.0

Rscript $script $geneTPM $output


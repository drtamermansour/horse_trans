ann="$1"
output="$2"
script="$3"

module load R/3.0.1

Rscript $script $ann $output

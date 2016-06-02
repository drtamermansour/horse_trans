ann="$1"
chromSizes="$2"
output="$3"
script="$4"

module load R/3.0.1

Rscript $script $ann $chromSizes $output

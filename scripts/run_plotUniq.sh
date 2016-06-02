summary0="$1"
summary5="$2"
output="$3"
expFolder="$4"
script="$5"

module load R/3.2.0

Rscript $script $summary0 $summary5 $output $expFolder


sup="$1"
cons="$2"
ORF="$3"
output="$4"
script="$5"

module load R/3.2.0

Rscript $script $sup $cons $ORF $output


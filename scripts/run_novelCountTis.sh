sup="$1"
cons="$2"
ORF="$3"
isoTPM="$4"
output="$5"
script="$6"

module load R/3.2.0

Rscript $script $sup $cons $ORF $isoTPM $output


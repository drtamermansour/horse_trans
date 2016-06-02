geneTPM="$1"
outputMap="$2"
outputTable="$3"
annTable="$4"
script="$5"

module load R/3.2.0

Rscript $script $geneTPM $outputMap $outputTable $annTable


genome="$1"
output="$2"
isoformfrac="$3"
assemblies="$4"

module load cufflinks/2.2.1

log=$(echo $output | sed 's/\//_/g').log
cuffmerge -s $genome -o $output --num-threads 4 --min-isoform-fraction $isoformfrac $assemblies > $log 2>&1
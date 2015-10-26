genome="$1"
output="$2"
assemblies="$3"

module load cufflinks/2.2.1

log=$(echo $output | sed 's/\//_/g').log
cuffmerge -s $genome -o $output --num-threads 4 $assemblies > $log 2>&1
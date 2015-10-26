genome="$1"
output="$2"
assemblies="$3"
Genes_GTF_file="$4"

module load cufflinks/2.2.1

log=$(echo $output | sed 's/\//_/g').log
cuffmerge -s $genome -o $output --num-threads 4 --ref-gtf $Genes_GTF_file $assemblies > $log 2>&1
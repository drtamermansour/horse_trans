refPtn="$1"

module load BLAST+/2.2.30
makeblastdb -in $refPtn -dbtype prot
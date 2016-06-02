module load Bioperl/1.6.923

trans="$1"
output="$2"
script_path=$(dirname "${BASH_SOURCE[0]}")  

perl $script_path/seq_stats.pl $trans > $output


genome="$1"

script_path=$(dirname "${BASH_SOURCE[0]}")

$script_path/UCSC_kent_commands/axtChain -linearGap=medium -psl NCBItoUCSC_map.psl -faT ncbi_genome2.fa -faQ $genome NCBItoUCSC_map.chain











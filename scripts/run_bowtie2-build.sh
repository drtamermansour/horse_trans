if [ $# -lt 3 ]
then
    printf "\nUsage run_bowtie2-build.sh [Genome file in fasta format] [Genome index base] [platform]\n"
    exit 0
fi

genome_file="$1"
genome_index_base="$2"
platform="$3"
script_path=$(dirname "${BASH_SOURCE[0]}")


if [ $platform == $"HPC" ];then
  bash ${script_path}/bowtie2-build.sh $genome_file $genome_index_base $platform
elif [ $platform == $"AMZ" ];then
  bash $script_path/bowtie2-build_AMZ.sh $genome_file $genome_index_base $platform
else echo "Platform should be either HPC or AMZ";fi



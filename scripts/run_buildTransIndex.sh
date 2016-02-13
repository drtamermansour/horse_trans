refGTF_file="$1"
transcriptome_index="$2"
Bowtie2_genome_index_base="$3"
platform="$4"                                                                                                                                            
script_path=$(dirname "${BASH_SOURCE[0]}")                                                                                                               

if [ $platform == $"HPC" ];then                                                                            
  bash ${script_path}/buildTransIndex.sh $refGTF_file $transcriptome_index $Bowtie2_genome_index_base
elif [ $platform == $"AMZ" ];then                                                                                                                      
  bash $script_path/buildTransIndex_AMZ.sh $refGTF_file $transcriptome_index $Bowtie2_genome_index_base                                              
else echo "Platform should be either HPC or AMZ";fi

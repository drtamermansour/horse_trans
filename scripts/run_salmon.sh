#!/bin/sh

if [ $# -lt 4 ]
then
printf "\nUsage run_tophat.sh [library type] [first strand] [transcriptome_index] [script_path]\n"
exit 0
fi

lib="$1"
strand="$2"
index="$3"
identifier="$4"
seq_dir="$5"
script_path="$6"

## determine lib type
## http://salmon.readthedocs.org/en/latest/library_type.html#fraglibtype
if [ $lib = $"PE" ];then orien="I"; script=$script_path/salmonQuant_PE.sh;
elif [ $lib = $"SE" ];then orien=""; script=$script_path/salmonQuant_SE.sh;
else exit; fi

if [ $strand = $"fr-unstranded" ];then st="U";
elif [ $strand = $"fr-firststrand" ];then st="SR";
elif [ $strand = $"fr-secondstrand" ];then st="SF";
else exit; fi

lib_type=$orien$st

qsub -v index="${index}",\
lib_type="${lib_type}",\
seq_dir="${seq_dir}",\
identifier="${identifier}" "$script"

if [ $lib = $"PE" ];then
qsub -v index="${index}",\
lib_type="U",\
seq_dir="${seq_dir}",\
identifier="${identifier}.SE" "$script_path/salmonQuant_SE.sh"
fi

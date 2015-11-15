#!/bin/sh

if [ $# -lt 3 ]
then
printf "\nUsage run_RealignerTargetCreator_forKnowns.sh [known indels] [indexed reference fasta] [script]\n"
exit 0
fi

indels="$1"
gatk_ref="$2"
script="$3"


qsub -v indels="${indels}",gatk_ref="${gatk_ref}" "${script}"

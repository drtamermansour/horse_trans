#!/bin/bash -login
#PBS -l walltime=02:00:00,nodes=1:ppn=1,mem=48Gb
#mdiag -A ged
#PBS -m abe
#PBS -N genome_to_cdna_fasta

#module load TransDecoder/2.0.1
#decoderUtil=$"/opt/software/TransDecoder/2.0.1--GCC-4.4.5/util"

cd $PBS_O_WORKDIR

##Construct the transcript fasta file
$decoderUtil/cufflinks_gtf_genome_to_cdna_fasta.pl $inputGTF $genome > $outputFASTA


qstat -f ${PBS_JOBID}

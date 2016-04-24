#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=4,mem=24Gb
#mdiag -A ged
#PBS -m abe
#PBS -N blastp

module load BLAST+/2.2.30

cd $PBS_O_WORKDIR

blastx -query $seq  -db $refPtn -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 4 > $output

qstat -f ${PBS_JOBID}

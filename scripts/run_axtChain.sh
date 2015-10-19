#!/bin/bash -login
#PBS -l walltime=1:00:00,nodes=1:ppn=1,mem=32Gb
#mdiag -A ged
#PBS -m abe
#PBS -N axtchain


cd $PBS_O_WORKDIR

$HOME/bin/UCSC_kent_commands/axtChain -linearGap=medium -psl NCBItoUCSC_map.psl -faT ncbi_genome2.fa -faQ $genome NCBItoUCSC_map.chain

qstat -f ${PBS_JOBID}







#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=1,mem=24Gb
#mdiag -A ged
#PBS -m abe
#PBS -N transdecoderPredict

module load TransDecoder/2.0.1

cd $PBS_O_WORKDIR

## predict the likely coding regions
TransDecoder.Predict -t $inputFASTA --retain_pfam_hits $pfamHomo --retain_blastp_hits $blastpHomo

qstat -f ${PBS_JOBID}

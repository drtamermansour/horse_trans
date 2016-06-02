#!/bin/bash -login
#PBS -l walltime=01:00:00,nodes=1:ppn=1,mem=1Gb
#mdiag -A ged
#PBS -m abe		#send email to myself
#PBS -N T_FastQC	#give name to the job

module load FastQC

cd $PBS_O_WORKDIR

fastqc -f fastq -noextract ${INPUT_FILE}

qstat -f ${PBS_JOBID}



## on the command line go to the data directory and run a loop for all qz files in this directory
# for f in *.gz; do qsub -v INPUT_FILE="$f" /mnt/home/mansourt/mRNAseq/fastqc.sh; done

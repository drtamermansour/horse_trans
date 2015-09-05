#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=1,mem=2Gb
#mdiag -A ged	
#PBS -m abe			#send email to myself
#PBS -N interleave		#give name to the job


cd $PBS_O_WORKDIR
source $HOME/env/bin/activate
module load GNU/4.7.1

interleave-reads.py $f1 $f2 > $output

qstat -f ${PBS_JOBID}


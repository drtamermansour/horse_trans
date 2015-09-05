#!/bin/bash -login
#mdiag -A ged
#PBS -m abe			#send email to myself
#PBS -N digi_norm		#give name to the job


source $HOME/env/bin/activate
module load GNU/4.7.1
cd $PBS_O_WORKDIR

#normalize-by-median.py -p -u $singletons -k $kmer -C $cutoff -N 4 -x $x --savegraph $kmer_table $samples
normalize-by-median.py -p -u $singletons -k $kmer -C $cutoff -N 4 -x $x $samples

qstat -f ${PBS_JOBID}


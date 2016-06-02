#!/bin/bash -login
#PBS -l walltime=2:00:00:00,nodes=1:ppn=16,mem=48Gb
#mdiag -A ged
#PBS -m abe
#PBS -N TransDecoder

module load TransDecoder/2.0.1
module load BLAST+/2.2.30
module load HMMER/3.0
decoderUtil=$"/opt/software/TransDecoder/2.0.1--GCC-4.4.5/util"

cd $PBS_O_WORKDIR

##Construct the transcript fasta file
$decoderUtil/cufflinks_gtf_genome_to_cdna_fasta.pl ../merged.gtf $genome > transcripts.fasta

## convert the transcript structure GTF file to an alignment-GFF3 formatted file
$decoderUtil/cufflinks_gtf_to_alignment_gff3.pl ../merged.gtf > transcripts.gff3

## extract the long open reading frames
TransDecoder.LongOrfs -t transcripts.fasta

## identify ORFs with homology to known proteins via blast and pfam searches.
blastp -query transcripts.fasta.transdecoder_dir/longest_orfs.pep  -db $refPtn  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 16 > blastp.outfmt6
hmmscan --cpu 16 --domtblout pfam.domtblout $refPfam transcripts.fasta.transdecoder_dir/longest_orfs.pep

## predict the likely coding regions
TransDecoder.Predict -t transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6

## generate a genome-based coding region annotation file
$decoderUtil/cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta 1> transcripts.fasta.transdecoder.genome.gff3 2> sterr
grep "Warning" sterr > warnings.txt && rm sterr

# convert the genome-based gene-gff3 file to bed
$decoderUtil/gff3_file_to_bed.pl transcripts.fasta.transdecoder.genome.gff3 > transcripts.fasta.transdecoder.genome.bed



qstat -f ${PBS_JOBID}

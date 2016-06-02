label="$1"

module load Trinotate/2.0.1     ## automatocally load trinity/20140413p1
## initiate the Trinotate database
Trinotate $label.Trinotate.sqlite init --gene_trans_map $label.gene_trans_map_file --transcript_fasta $label.fasta --transdecoder_pep $label.pep
## Loading BLAST homologies
Trinotate $label.Trinotate.sqlite LOAD_swissprot_blastp $label.blastp.outfmt6
Trinotate $label.Trinotate.sqlite LOAD_swissprot_blastx $label.blastx.outfmt6
##  Load Pfam domain entries
Trinotate $label.Trinotate.sqlite LOAD_pfam $label.pfam.domtblout
## Load transmembrane domains
#Trinotate $label.Trinotate.sqlite LOAD_tmhmm $label.tmhmm.out
## Load signal peptide predictions
#Trinotate $label.Trinotate.sqlite LOAD_signalp $label.signalp.out
## Load RNAMMER predictions
#Trinotate $label.Trinotate.sqlite LOAD_rnammer $label.rnammer.gff
## Output an Annotation Report
Trinotate $label.Trinotate.sqlite report > $label.annotation_report.xls


myRoot=$"/mnt/ls12/Tamer"
resources=$myRoot/horse_trans/resources
##################################################################################################

mkdir -p $resources/NCBI/GFF
cd $resources/NCBI/GFF
wget ftp://ftp.ncbi.nih.gov/genomes/Equus_caballus/GFF/ref_EquCab2.0_top_level.gff3.gz
gunzip ref_EquCab2.0_top_level.gff3.gz
cat ref_EquCab2.0_top_level.gff3 | awk  -F $'\t' '{print $2,"\t",$3}' | sort | uniq > id_type
grep -v "^#" ref_EquCab2.0_top_level.gff3 > ref_EquCab2.0_top_level.gff3.nohash
cat ref_EquCab2.0_top_level.gff3.nohash | awk -F $'\t' '{A[$2"\t"$3]++}END{for(i in A)print i,A[i]}' | sort > id_type.count

grep -v "NC_001640.1" ref_EquCab2.0_top_level.gff3.nohash > ref_EquCab2.0_top_level.gff3.nohash.noMT
grep "NC_001640.1" ref_EquCab2.0_top_level.gff3.nohash > ref_EquCab2.0_top_level.gff3.nohash.MT

cat ref_EquCab2.0_top_level.gff3.nohash.noMT | awk -F $'\t' '{A[$2"\t"$3]++}END{for(i in A)print i,A[i]}' | sort > id_type.count2
cat ref_EquCab2.0_top_level.gff3.nohash.MT | awk -F $'\t' '{A[$2"\t"$3]++}END{for(i in A)print i,A[i]}' | sort > id_type.count3

grep "gbkey=Gene;" ref_EquCab2.0_top_level.gff3 | wc -l
cat ref_EquCab2.0_top_level.gff3 | awk -F $'\t' '$3 == "mRNA"'| wc -l
cat ref_EquCab2.0_top_level.gff3 | awk -F $'\t' '$3 == "mRNA"' | awk -F $';' '{print $3}' | sort | uniq | wc -l
## for the new annotation 11/25/2015
#cat ref_EquCab2.0_top_level.gff3 | awk -F $'\t' '$3 == "mRNA"' | awk -F $';' '{print $2}' | sort | uniq | wc -l

## pseudo genes
grep pseudo=true ref_EquCab2.0_top_level.gff3 | wc -l

cat ref_EquCab2.0_top_level.gff3.nohash | awk -F $'\t' '$3 == "cDNA_match"' > cDNA_match
cat cDNA_match | awk -F $'\t' '{print $9}' | awk -F $';' '{for(i=1;i<=NF;i++)printf "%s ",$i; print ""}' | sort -u -k1,1 --merge > uniq_cDNA_match
wc -l uniq_cDNA_match
grep "Target=NM_" uniq_cDNA_match | wc -l
grep "Target=XM_" uniq_cDNA_match | wc -l
cat ref_EquCab2.0_top_level.gff3.nohash | awk -F $'\t' '$3 == "match"' > match
cat match | awk -F $'\t' '{print $9}' | awk -F $';' '{for(i=1;i<=NF;i++)printf "%s ",$i; print ""}' | sort -u -k1,1 --merge > uniq_match
wc -l uniq_match

mkdir -p $resources/NCBI/RNA
cd $resources/NCBI/RNA
wget ftp://ftp.ncbi.nih.gov/genomes/Equus_caballus/RNA/rna.fa.gz
gunzip rna.fa.gz
##################################################################################################

mkdir -p $resources/UCSC_tables
cd $resources/UCSC_tables
cat RefSeq.GTF | awk -F $'\t' '{print $9}' | awk '{print $4}' | awk -F $'\"' '{print $2}' | sort | uniq > transcripts
wc -l transcripts   ## 1889
grep "^NM_" transcripts > transcripts_NM
grep "^NR_" transcripts > transcripts_NR
cat transcripts_NM | awk '{print $1}' | while read pat; do grep "^$pat" refseq_quary_NM; done > matches
## create the refseq_quary files from quaries of Entrez website with selection for mRNA and other abaliable RNAs
cat refseq_quary_NM | awk -F $'.' '{print $1}' > refseq_quary.Acc_NM
cat refseq_quary_NR | awk -F $'.' '{print $1}' > refseq_quary.Acc_NR
cat refseq_quary_XM | awk -F $'.' '{print $1}' > refseq_quary.Acc_XM
## create the refseq_quary files from quaries of Entrez website without selection
grep "^NM" refseq_quary_ALL > refseq_quary_NMext
grep "^XM" refseq_quary_ALL > refseq_quary_XMext
grep "^NR" refseq_quary_ALL > refseq_quary_NRext
grep "^XR" refseq_quary_ALL > refseq_quary_XRext
grep "^NG" refseq_quary_ALL > refseq_quary_NGext
grep "^NC" refseq_quary_ALL > refseq_quary_genome
grep "^NW" refseq_quary_ALL >> refseq_quary_genome
grep "^NT" refseq_quary_ALL >> refseq_quary_genome
ls refseq_quary_*ext | xargs  wc -l
cat refseq_quary_NMext | awk -F $'.' '{print $1}' > refseq_quary.Acc_NMext
cat refseq_quary_XMext | awk -F $'.' '{print $1}' > refseq_quary.Acc_XMext
cat refseq_quary_NRext | awk -F $'.' '{print $1}' > refseq_quary.Acc_NRext
cat refseq_quary_XRext | awk -F $'.' '{print $1}' > refseq_quary.Acc_XRext

## grep -F -x -v -f fileB fileA  ## printing out the lines in fileA that don't contain the same data as any line in fileB.
grep -F -x -v -f transcripts_NM refseq_quary.Acc_NMext > refseq_quary.Acc_NMext.uniq
grep -F -x -v -f transcripts_NR refseq_quary.Acc_NRext > refseq_quary.Acc_NRext.uniq
## OR ## comm <(sort transcripts_NM) <(sort refseq_quary.Acc_NM) -3 > refseq_quary.Acc_NM.uniq
grep -F -x -v -f refseq_quary.Acc_NMext transcripts_NM > transcripts_NM.uniq
grep -F -x -v -f refseq_quary.Acc_NRext transcripts_NR > transcripts_NR.uniq
#########################################################################################################

mkdir -p $resources/ensemble
cd $resources/ensemble
wget ftp://ftp.ensembl.org/pub/release-80/gtf/equus_caballus/*
gunzip *.gz
grep -v "^#" Equus_caballus.EquCab2.80.gtf > Equus_caballus.EquCab2.80.gtf.nohashq
cat Equus_caballus.EquCab2.80.gtf.nohash | awk -F $'\t' '{A[$3]++}END{for(i in A)print i,A[i]}' | sort > id_type.count
cat Equus_caballus.EquCab2.80.gtf.nohash | awk -F $'\t' '$3 == "transcript"' | awk -F $'\t' '{print $9}' | awk '{A[$NF]++}END{for(i in A)print i,A[i]}' | sort > gene_types.count
# Using Biomart website, create the list of ensemble genes and corresponding RefSeq IDs
grep -v "^Ensembl.Gene.ID" ensamble_biomart.tab > ensamble_biomart.tab.noheader
cat ensamble_biomart.tab.noheader | awk -F $'\t' '{print $13}' | sed '/^\s*$/d' | wc -l
cat ensamble_biomart.tab.noheader | awk -F $'\t' '{print $13}' | sed '/^\s*$/d' | sort | uniq | wc -l
cat ensamble_biomart.tab.noheader | awk -F $'\t' '{print $14}' | sed '/^\s*$/d' | wc -l
cat ensamble_biomart.tab.noheader | awk -F $'\t' '{print $14}' | sed '/^\s*$/d' | sort | uniq | wc -l
cat ensamble_biomart.tab.noheader | awk -F $'\t' '{print $15}' | sed '/^\s*$/d' | wc -l
cat ensamble_biomart.tab.noheader | awk -F $'\t' '{print $15}' | sed '/^\s*$/d' | sort | uniq | wc -l
cat ensamble_biomart.tab.noheader | awk -F $'\t' '{print $16}' | sed '/^\s*$/d' | wc -l
cat ensamble_biomart.tab.noheader | awk -F $'\t' '{print $16}' | sed '/^\s*$/d' | sort | uniq | wc -l
cat ensamble_biomart.tab.noheader | awk -F $'\t' '$13 == ""' | awk -F $'\t' '{print $15}' | sed '/^\s*$/d' | sort | uniq | wc -l
cat ensamble_biomart.tab.noheader | awk -F $'\t' '$15 == ""' | awk -F $'\t' '{print $13}' | sed '/^\s*$/d' | sort | uniq | wc -l
cat ensamble_biomart.tab.noheader | awk -F $'\t' '$13 != ""' | awk -F $'\t' '$15 != ""' | awk -F $'\t' '{print $13'\t'$15}' | sed '/^\s*$/d' | sort | uniq | wc -l
cat ensamble_biomart.tab.noheader | awk -F $'\t' '($13 == "" || $14 == "" || $15 == "" || $16 == "")' | awk -F $'\t' '{print $2}' | sed '/^\s*$/d' | sort | uniq | wc -l

cat ensamble_biomart.tab.noheader | awk -F $'\t' '($13 != "" || $14 != "" || $15 != "" || $16 != "")' | awk -F $'\t' '$11 == "NOVEL"' | wc -l


cat ensamble_biomart.tab.noheader | awk -F $'\t' '$4 == "protein_coding"' > ensamble_biomart.tab.noheader.ptn
cat ensamble_biomart.tab.noheader.ptn | awk -F $'\t' '{print $1}' | sort | uniq | wc -l ## no of protein coding genes






















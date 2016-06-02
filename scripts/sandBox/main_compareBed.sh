## compare bed files
mkdir $horse_trans/compareBed
cd $horse_trans/compareBed

## define the list of home made annoatations
alltissueGTF_file=$tissue_Cuffmerge/$"all_tissues"/$cuffmerge_run/merged.gtf

## get a copy of the annotation you want to compare
for f in refGTF_file ncbiNoNameGTF_file ensGTF_file alltissueGTF_file; do
  cp ${!f} ${f}.gtf
  cat ${!f} | gzip > ${f}.gtf.gz
done
## create ExonMerge trackhup
## merge transcripts per loci
#for f in *GTF_file.gtf; do
#  filename=$(basename ${f%.gtf})
#  cat $f \
#    | cgat gtf2gtf --method=sort --sort-order=gene \
#    | cgat gtf2gtf --method=merge-exons --with-utr \
#    | cgat gtf2gtf --method=set-transcript-to-gene \
#    | cgat gtf2gtf --method=sort --sort-order=position \
#    > ${filename}_mergeExons_withUTR.gtf
#done
for f in ncbiNoNameGTF_file ensGTF_file alltissueGTF_file; do
  filename=$(basename ${f%.gtf})
  cat $f \
    | cgat gtf2gtf --method=sort --sort-order=gene \
    | cgat gtf2gtf --method=merge-exons \
    | cgat gtf2gtf --method=set-transcript-to-gene \
    | cgat gtf2gtf --method=sort --sort-order=position \
    > ${filename}_mergeExons.gtf
done
## create list of ExonMerge assemblies
## & convert reference GTF to bed files
rm -f exonMerge_assemblies.txt
for f in *_mergeExons.gtf; do
  bash $script_path/gtfToBigBed.sh "$f" "$genome_dir/$UCSCgenome.chrom.sizes" "$script_path"
  identifier=${f%.gtf}
  cp ${identifier}.BigBed $track_hub/$UCSCgenome/BigBed/.
  echo ${identifier} >> exonMerge_assemblies.txt;
done
## initiate a given track hub
hub_name=$"HorseTrans_exonMerge_assemblies"
shortlabel=$"exonMerge_assemblies"
longlabel=$"Assemblies presented with exon merge"
email=$"drtamermansour@gmail.com"
cd $track_hub
bash $script_path/create_trackHub.sh "$UCSCgenome" "$hub_name" "$shortlabel" "$longlabel" "$email"
## edit the trackDb
> $horse_trans/emptyTemp.txt
current_libs=$track_hub/current_libs_$shortlabel
current_tissues=$track_hub/current_tiss_$shortlabel
trackDb=$track_hub/$UCSCgenome/trackDb_$shortlabel.txt
lib_assemblies=$horse_trans/compareBed/exonMerge_assemblies.txt
tiss_assemblies=$horse_trans/emptyTemp.txt
bash $script_path/edit_trackDb.sh $current_libs $current_tissues $trackDb $lib_assemblies $tiss_assemblies
## compare gene sets
python /CAGT_devel/cgat-code/scripts/diff_gtf.py alltissueGTF_file.gtf ncbiNoNameGTF_file.gtf ensGTF_file.gtf > threeway.tsv
python /CAGT_devel/cgat-code/scripts/diff_gtf.py --update=threeway.tsv  alltissueGTF_file.gtf ncbiNoNameGTF_file.gtf ensGTF_file.gtf refGTF_file.gtf > fourway.tsv

mkdir allVSncbi && cd allVSncbi
python /CAGT_devel/cgat-code/scripts/gtfs2tsv.py ../alltissueGTF_file.gtf ../ncbiNoNameGTF_file.gtf > allVSncbi.tsv
mkdir ../allVSens && cd ../allVSens
python /CAGT_devel/cgat-code/scripts/gtfs2tsv.py ../alltissueGTF_file.gtf ../ensGTF_file.gtf > allVSens.tsv
mkdir ../ensVSncbi && cd ../ensVSncbi
python /CAGT_devel/cgat-code/scripts/gtfs2tsv.py ../ensGTF_file.gtf ../ncbiNoNameGTF_file.gtf > ensVSncbi.tsv

## do the bed comparison
filename=${alltissueGTF_file%.gtf}
$script_path/UCSC_kent_commands/gtfToGenePred $alltissueGTF_file ${filename}.gpred
cat ${filename}.gpred | $script_path/genePredToBed > ${filename}.bed
alltissueBED_file=${filename}.bed

for f in alltissue ncbiNoName ens ref;do
echo $f >> assAnn.txt
bed="$f"BED_file
echo "No_of_transcripts=" $(wc -l ${!bed}) >> assAnn.txt       ## No of transcripts 157567
cat ${!bed} | awk -F $'\t' '{A["Transcripts with "$10" exons= "]++}END{for(i in A)print i,A[i]}' | sort -n > $f.exonsPerTranscript.count
echo "Trans_W_1_exon=" $(cat ${!bed} | awk -F $'\t' '$10 == 1' | wc -l) >> assAnn.txt ## Trans w 1 exon=  58571
echo "Trans_W_2_exon=" $(cat ${!bed} | awk -F $'\t' '$10 == 2' | wc -l) >> assAnn.txt ## Trans w 2 exons=  8624
echo "Trans_W_>2_exon=" $(cat ${!bed} | awk -F $'\t' '$10 > 1' | wc -l) >> assAnn.txt ## Trans w >2 exons=90372

mergeExon="$f"GTF_file_mergeExons.bed
echo "No_of_genes=" $(wc -l $mergeExon) >> assAnn.txt     ## No of genes 73808
cat $mergeExon | awk -F $'\t' '{A["Gene with "$10" exons= "]++}END{for(i in A)print i,A[i]}' | sort -n > $f.exonsPerGene.count
echo "Genes_W_1_exon=" $(cat $mergeExon | awk -F $'\t' '$10 == 1' | wc -l) >> assAnn.txt ## genes w 1 exon=54822
echo "Genes_W_2_exon=" $(cat $mergeExon | awk -F $'\t' '$10 == 1' | wc -l) >> assAnn.txt ## genes w 2 exons=2707
echo "Genes_W_>2_exon=" $(cat $mergeExon | awk -F $'\t' '$10 == 1' | wc -l) >> assAnn.txt ## genes w >2 exons=90372
done


#cat $alltissueBED_file | awk -F $'\t' '{A[$2"\t"$3]++}END{for(i in A)print i,A[i]}' | sort > id_type.count
module load cufflinks/2.2.1
gffread -E $alltissueGTF_file -o ${filename}.gff
alltissueGFF_file=${filename}.gff

module load BEDTools/2.24.0


script_path=$myRoot/horse_trans/scripts
rawData=$myRoot/horse_trans/rawdata
prepData=$myRoot/horse_trans/prepdata

###############################################################################################
mkdir -p $rawData/Muscle/PE_125_fr.firststrand_Valberg.Finno_01012015/rawdata
cd $rawData/Muscle/PE_125_fr.firststrand_Valberg.Finno_01012015/rawdata
rsync -avzP -e ssh "tmansour@kaiden.genomecenter.ucdavis.edu:/finno/data/galaxy_data/WB_RNASeq/WB_RNASeq/raw_data/*" .
bash ${script_path}/run_fastqc.sh "$rawData/Muscle/PE_125_fr.firststrand_Valberg.Finno_01012015/rawdata"
mkdir -p $prepData/Muscle/PE_125_fr.firststrand_Valberg.Finno_01012015/fastq_data
mv $rawData/Muscle/PE_125_fr.firststrand_Valberg.Finno_01012015/rawdata/* $prepData/Muscle/PE_125_fr.firststrand_Valberg.Finno_01012015/fastq_data/.
###############################################################################################
mkdir -p $rawData/BrainStem/PE_100_fr.firststrand_Finno_01012015/rawdata
cd $rawData/BrainStem/PE_100_fr.firststrand_Finno_01012015/rawdata
rsync -avzP -e ssh "tmansour@kaiden.genomecenter.ucdavis.edu:/finno/data/galaxy_data/NAD_RNASeq_Bstem_July2013/raw_reads/*.fastq.gz" .
bash ${script_path}/run_fastqc.sh "$rawData/BrainStem/PE_100_fr.firststrand_Finno_01012015/rawdata"
mkdir -p $prepData/BrainStem/PE_100_fr.firststrand_Finno_01012015/fastq_data
mv $rawData/BrainStem/PE_100_fr.firststrand_Finno_01012015/rawdata/* $prepData/BrainStem/PE_100_fr.firststrand_Finno_01012015/fastq_data/.
###############################################################################################
mkdir -p $rawData/SpinalCord/PE_100_fr.firststrand_Finno_01012015/rawdata
cd $rawData/SpinalCord/PE_100_fr.firststrand_Finno_01012015/rawdata
rsync -avzP -e ssh "tmansour@kaiden.genomecenter.ucdavis.edu:/finno/data/galaxy_data/NAD_RNASeq_SC_Dec2013/fastq/raw_data/*_SC*.fastq.gz" .
bash ${script_path}/run_fastqc.sh "$rawData/SpinalCord/PE_100_fr.firststrand_Finno_01012015/rawdata"
mkdir -p $prepData/SpinalCord/PE_100_fr.firststrand_Finno_01012015/fastq_data
mv $rawData/SpinalCord/PE_100_fr.firststrand_Finno_01012015/rawdata/* $prepData/SpinalCord/PE_100_fr.firststrand_Finno_01012015/fastq_data/.
###############################################################################################
mkdir -p $rawData/Cerebellum/PE_100_fr.firststrand_Scott.Murray.Penedo_01012015/rawdata
cd $rawData/Cerebellum/PE_100_fr.firststrand_Scott.Murray.Penedo_01012015/rawdata
rsync -avzP -e ssh "tmansour@kaiden.genomecenter.ucdavis.edu:/finno/data/eyscottRNAseq/data/RawReads/*.fastq.gz" .
rsync -avzP -e ssh "tmansour@kaiden.genomecenter.ucdavis.edu:/finno/data/eyscottRNAseq/data/RawReads/2-3_MO_C_GATCAG_L007_R1_001.fastq.zip" .
unzip 2-3_MO_C_GATCAG_L007_R1_001.fastq.zip
rm -r __MACOSX
rm 2-3_MO_C_GATCAG_L007_R1_001.fastq.zip
gzip 2-3_MO_C_GATCAG_L007_R1_001.fastq
> replicates.txt
for f in *_L006_R1_001.fastq.gz; do f2=$(basename $f | sed 's/_L006_R1_/_L007_R1_/'); echo $f $f2 >> replicates.txt; done
bash ${script_path}/run_fastqc.sh "$rawData/Cerebellum/PE_100_fr.firststrand_Scott.Murray.Penedo_01012015/rawdata"
mkdir -p $prepData/Cerebellum/PE_100_fr.firststrand_Scott.Murray.Penedo_01012015/fastq_data
mv $rawData/Cerebellum/PE_100_fr.firststrand_Scott.Murray.Penedo_01012015/rawdata/* $prepData/Cerebellum/PE_100_fr.firststrand_Scott.Murray.Penedo_01012015/fastq_data/.
###############################################################################################
mkdir -p $rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/rawdata
cd $rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/rawdata
rsync -avzP -e ssh "tmansour@kaiden.genomecenter.ucdavis.edu:/finno/data/Bellone_RNAseq/reads/{10220687_61G3EAAXX_s_6_1_sequence.txt,10220687_61G3EAAXX_s_6_2_sequence.txt}" .
module load FASTX/0.0.14
mkdir -p $rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp1
cat 10220687_61G3EAAXX_s_6_1_sequence.txt | fastx_barcode_splitter.pl --bcfile barcodes_for_7.54_7.45.txt --bol --mismatches 1 --prefix $rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp1/Retina_ --suffix "_R1_001.fastq"
cat 10220687_61G3EAAXX_s_6_2_sequence.txt | fastx_barcode_splitter.pl --bcfile barcodes_for_7.54_7.45.txt --bol --mismatches 1 --prefix $rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp1/Retina_ --suffix "_R2_001.fastq"
rm $rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp1/*unmatched*.fastq

cd $rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp1
mkdir -p $rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp2
for f in *.fastq; do echo $f; fastx_trimmer -f 6 -i $f -o $rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp2/$f; done
## fix the broken paired ends
cd $rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp2
mkdir -p $rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp3
for f in *_R1_*.fastq; do
    f2=$(basename $f | sed 's/_R1_/_R2_/')
    fs=$(basename $f | sed 's/_R1_/_R1S_/')
    f2s=$(basename $f | sed 's/_R1_/_R2S_/')
    $HOME/perl5/bin/pairfq makepairs -f $f -r $f2 \
-fp $rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp3/$f \
-rp $rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp3/$f2 \
-fs $rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp3/$fs \
-rs $rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp3/$f2s
done
cd $rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp3
for f in *.fastq; do echo $f; gzip $f; done
bash ${script_path}/run_fastqc.sh "$rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp3"
#QC check shows that the data has encoding illumina 1.5 so we need to change this to sanger encoding
gunzip *.fastq.gz
mkdir -p $rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp4
module load EMBOSS/6.5.7
for f in *_R[1-2]_*.fastq; do seqret fastq-illumina::$f fastq::$rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp4/$f; done
cd $rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp4
for f in *.fastq; do gzip $f; done
bash ${script_path}/run_fastqc.sh "$rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp4"
mkdir -p $prepData/Retina/PE_81_fr.unstranded_Bellone_01012015/fastq_data
mv $rawData/Retina/PE_81_fr.unstranded_Bellone_01012015/temp4/* $prepData/Retina/PE_81_fr.unstranded_Bellone_01012015/fastq_data/.
###############################################################################################
mkdir -p $rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/rawdata
cd $rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/rawdata
rsync -avzP -e ssh "tmansour@kaiden.genomecenter.ucdavis.edu:/finno/data/Bellone_RNAseq/reads/{10220687_61G3EAAXX_s_7_1_sequence.txt,10220687_61G3EAAXX_s_7_2_sequence.txt}" .
module load FASTX/0.0.14
mkdir -p $rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp1
cat 10220687_61G3EAAXX_s_7_1_sequence.txt | fastx_barcode_splitter.pl --bcfile barcodes_for_7.49.txt --bol --mismatches 1 --prefix $rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp1/Skin_ --suffix "_R1_001.fastq"
cat 10220687_61G3EAAXX_s_7_2_sequence.txt | fastx_barcode_splitter.pl --bcfile barcodes_for_7.49.txt --bol --mismatches 1 --prefix $rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp1/Skin_ --suffix "_R2_001.fastq"
rm $rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp1/*unmatched*.fastq

cd $rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp1
mkdir -p $rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp2
for f in *.fastq; do echo $f; fastx_trimmer -f 6 -i $f -o $rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp2/$f; done
## fix the broken paired ends
cd $rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp2
mkdir -p $rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp3
for f in *_R1_*.fastq; do
    f2=$(basename $f | sed 's/_R1_/_R2_/')
    fs=$(basename $f | sed 's/_R1_/_R1S_/')
    f2s=$(basename $f | sed 's/_R1_/_R2S_/')
    $HOME/perl5/bin/pairfq makepairs -f $f -r $f2 \
-fp $rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp3/$f \
-rp $rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp3/$f2 \
-fs $rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp3/$fs \
-rs $rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp3/$f2s
done
cd $rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp3
for f in *.fastq; do echo $f; gzip $f; done
bash ${script_path}/run_fastqc.sh "$rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp3"
#QC check shows that the data has encoding illumina 1.5 so we need to change this to sanger encoding
gunzip *.fastq.gz
mkdir -p $rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp4
module load EMBOSS/6.5.7
for f in *_R[1-2]_*.fastq; do seqret fastq-illumina::$f fastq::$rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp4/$f; done
cd $rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp4
for f in *.fastq; do gzip $f; done
bash ${script_path}/run_fastqc.sh "$rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp4"
mkdir -p $prepData/Skin/PE_81_fr.unstranded_Bellone_01012015/fastq_data
mv $rawData/Skin/PE_81_fr.unstranded_Bellone_01012015/temp4/* $prepData/Skin/PE_81_fr.unstranded_Bellone_01012015/fastq_data/.
###############################################################################################
mkdir -p $rawData/Skin/SE_81_fr.unstranded_Bellone_01012015/rawdata
cd $rawData/Skin/SE_81_fr.unstranded_Bellone_01012015/rawdata
rsync -avzP -e ssh "tmansour@kaiden.genomecenter.ucdavis.edu:/finno/data/Bellone_RNAseq/reads/10216429_61GJ3AAXX_s_6_sequence.txt" .
module load FASTX/0.0.14
mkdir -p $rawData/Skin/SE_81_fr.unstranded_Bellone_01012015/temp1
cat 10216429_61GJ3AAXX_s_6_sequence.txt | fastx_barcode_splitter.pl --bcfile barcodes_for_7.54_7.46.txt --bol --mismatches 1 --prefix $rawData/Skin/SE_81_fr.unstranded_Bellone_01012015/temp1/Skin_ --suffix "_SR_001.fastq"
cat 10216429_61GJ3AAXX_s_6_sequence.txt | fastx_barcode_splitter.pl --bcfile barcodes_for_7.54_7.46.txt --bol --mismatches 1 --prefix $rawData/Skin/SE_81_fr.unstranded_Bellone_01012015/temp1/Skin_ --suffix "_SR_001.fastq"
rm $rawData/Skin/SE_81_fr.unstranded_Bellone_01012015/temp1/*unmatched*.fastq

cd $rawData/Skin/SE_81_fr.unstranded_Bellone_01012015/temp1
mkdir -p $rawData/Skin/SE_81_fr.unstranded_Bellone_01012015/temp2
for f in *.fastq; do echo $f; fastx_trimmer -f 6 -i $f -o $rawData/Skin/SE_81_fr.unstranded_Bellone_01012015/temp2/$f; done
cd $rawData/Skin/SE_81_fr.unstranded_Bellone_01012015/temp2
for f in *.fastq; do echo $f; gzip $f; done
bash ${script_path}/run_fastqc.sh "$rawData/Skin/SE_81_fr.unstranded_Bellone_01012015/temp2"
#QC check shows that the data has encoding illumina 1.5 so we need to change this to sanger encoding
gunzip *.fastq.gz
mkdir -p $rawData/Skin/SE_81_fr.unstranded_Bellone_01012015/temp3
module load EMBOSS/6.5.7
for f in *_SR_*.fastq; do seqret fastq-illumina::$f fastq::$rawData/Skin/SE_81_fr.unstranded_Bellone_01012015/temp3/$f; done
cd $rawData/Skin/SE_81_fr.unstranded_Bellone_01012015/temp3
for f in *.fastq; do gzip $f; done
bash ${script_path}/run_fastqc.sh "$rawData/Skin/SE_81_fr.unstranded_Bellone_01012015/temp3"
mkdir -p $prepData/Skin/SE_81_fr.unstranded_Bellone_01012015/fastq_data
mv $rawData/Skin/SE_81_fr.unstranded_Bellone_01012015/temp3/* $prepData/Skin/SE_81_fr.unstranded_Bellone_01012015/fastq_data/.

###############################################################################################
mkdir -p $rawData/Skin/SE_95_fr.unstranded_Bellone_01012015/rawdata
cd $rawData/Skin/SE_95_fr.unstranded_Bellone_01012015/rawdata
rsync -avzP -e ssh "tmansour@kaiden.genomecenter.ucdavis.edu:/finno/data/Bellone_RNAseq/reads/10237118_B08AEABXX_s_*" .
module load FASTX/0.0.14
mkdir -p $rawData/Skin/SE_95_fr.unstranded_Bellone_01012015/temp1
cp 10237118_B08AEABXX_s_7_D-052A_sequence.txt $rawData/Skin/SE_95_fr.unstranded_Bellone_01012015/temp1/Skin_Bellone_d052Pig_SR_001.fastq
cp 10237118_B08AEABXX_s_7_06-92A-pig_sequence.txt $rawData/Skin/SE_95_fr.unstranded_Bellone_01012015/temp1/Skin_Bellone_692Pig_SR_001.fastq
cp 10237118_B08AEABXX_s_7_06-92A-unpig_sequence.txt $rawData/Skin/SE_95_fr.unstranded_Bellone_01012015/temp1/Skin_Bellone_692NonPig_SR_001.fastq

cd $rawData/Skin/SE_95_fr.unstranded_Bellone_01012015/temp1
mkdir -p $rawData/Skin/SE_95_fr.unstranded_Bellone_01012015/temp2
for f in *.fastq; do echo $f; fastx_trimmer -f 6 -i $f -o $rawData/Skin/SE_95_fr.unstranded_Bellone_01012015/temp2/$f; done
cd $rawData/Skin/SE_95_fr.unstranded_Bellone_01012015/temp2
for f in *.fastq; do echo $f; gzip $f; done
bash ${script_path}/run_fastqc.sh "$rawData/Skin/SE_95_fr.unstranded_Bellone_01012015/temp2"
#QC check shows that the data has encoding illumina 1.5 so we need to change this to sanger encoding
gunzip *.fastq.gz
mkdir -p $rawData/Skin/SE_95_fr.unstranded_Bellone_01012015/temp3
module load EMBOSS/6.5.7
for f in *_SR_*.fastq; do seqret fastq-illumina::$f fastq::$rawData/Skin/SE_95_fr.unstranded_Bellone_01012015/temp3/$f; done
cd $rawData/Skin/SE_95_fr.unstranded_Bellone_01012015/temp3
for f in *.fastq; do gzip $f; done
bash ${script_path}/run_fastqc.sh "$rawData/Skin/SE_95_fr.unstranded_Bellone_01012015/temp3"
mkdir -p $prepData/Skin/SE_95_fr.unstranded_Bellone_01012015/fastq_data
mv $rawData/Skin/SE_95_fr.unstranded_Bellone_01012015/temp3/* $prepData/Skin/SE_95_fr.unstranded_Bellone_01012015/fastq_data/.
###############################################################################################
## Cell lineages of mammalian pre-implantation development: trophectoderm (TE), epiblast (EPI) & primitive endoderm (PE)
## The first cell fate decision segregates the TE from the inner cell mass (ICM).
## Prior to implantation, ICM gives rise to PE, which is a monolayer separating the blastocoel from the cluster of pluripotent EPI cells.
## The EPI forms the future fetus, the TE develops into the fetal placenta, and the PE becomes the visceral and parietal endoderm of the yolk sacs.
## The ICM/EPI is the source of pluripotent embryonic stem cells (ESCs)

mkdir -p $rawData/Embryo.ICM/SE_100_fr.unstranded_Ross_01032012/rawdata
cd $rawData/Embryo.ICM/SE_100_fr.unstranded_Ross_01032012/rawdata
rsync -avzP -e ssh "tmansour@kaiden.genomecenter.ucdavis.edu:/finno/data/Ross_Data/Equine\ ICM-TE\ RNAseq\ March\ 2012/{Sample_JEQ_index10,Sample_JEQ_index1,Sample_JEQ_index7}" .
mkdir -p $rawData/Embryo.ICM/SE_100_fr.unstranded_Ross_01032012/temp1
cat Sample_JEQ_index10/*_R1_*.fastq.gz > $rawData/Embryo.ICM/SE_100_fr.unstranded_Ross_01032012/temp1/ICM22_JEQ1_index10_SR_001.fastq.gz
cat Sample_JEQ_index1/*_R1_*.fastq.gz > $rawData/Embryo.ICM/SE_100_fr.unstranded_Ross_01032012/temp1/ICM15_JEQ3_index1_SR_001.fastq.gz
cat Sample_JEQ_index7/*_R1_*.fastq.gz > $rawData/Embryo.ICM/SE_100_fr.unstranded_Ross_01032012/temp1/ICM25_JEQ4_index7_SR_001.fastq.gz
bash ${script_path}/run_fastqc.sh "$rawData/Embryo.ICM/SE_100_fr.unstranded_Ross_01032012/temp1"
mkdir -p $prepData/Embryo.ICM/SE_100_fr.unstranded_Ross_01032012/fastq_data
mv $rawData/Embryo.ICM/SE_100_fr.unstranded_Ross_01032012/temp1/* $prepData/Embryo.ICM/SE_100_fr.unstranded_Ross_01032012/fastq_data/.

mkdir -p $rawData/Embryo.TE/SE_100_fr.unstranded_Ross_01032012/rawdata
cd $rawData/Embryo.TE/SE_100_fr.unstranded_Ross_01032012/rawdata
rsync -avzP -e ssh "tmansour@kaiden.genomecenter.ucdavis.edu:/finno/data/Ross_Data/Equine\ ICM-TE\ RNAseq\ March\ 2012/{Sample_JEQ_index11,Sample_JEQ_index8,Sample_JEQ_index9}" .
mkdir -p $rawData/Embryo.TE/SE_100_fr.unstranded_Ross_01032012/temp1
cat Sample_JEQ_index11/*_R1_*.fastq.gz > $rawData/Embryo.TE/SE_100_fr.unstranded_Ross_01032012/temp1/TE2412_JEQ2_index11_SR_001.fastq.gz
cat Sample_JEQ_index8/*_R1_*.fastq.gz > $rawData/Embryo.TE/SE_100_fr.unstranded_Ross_01032012/temp1/TE16_JEQ5_index8_SR_001.fastq.gz
cat Sample_JEQ_index9/*_R1_*.fastq.gz > $rawData/Embryo.TE/SE_100_fr.unstranded_Ross_01032012/temp1/TE2312_JEQ6_index9_SR_001.fastq.gz
bash ${script_path}/run_fastqc.sh "$rawData/Embryo.TE/SE_100_fr.unstranded_Ross_01032012/temp1"
mkdir -p $prepData/Embryo.TE/SE_100_fr.unstranded_Ross_01032012/fastq_data
mv $rawData/Embryo.TE/SE_100_fr.unstranded_Ross_01032012/temp1/* $prepData/Embryo.TE/SE_100_fr.unstranded_Ross_01032012/fastq_data/.

mkdir -p $rawData/Embryo.ICM/PE_100_fr.unstranded_Ross_01022013/rawdata
cd $rawData/Embryo.ICM/PE_100_fr.unstranded_Ross_01022013/rawdata
rsync -avzP -e ssh "tmansour@kaiden.genomecenter.ucdavis.edu:/finno/data/Ross_Data/Project_Ross_P/{Sample_PRKI001A_index1,Sample_PRKI001B_index3,Sample_PRKI001C_index4}" .
mkdir -p $rawData/Embryo.ICM/PE_100_fr.unstranded_Ross_01022013/temp1
cat Sample_PRKI001A_index1/*_R1_*.fastq.gz > $rawData/Embryo.ICM/PE_100_fr.unstranded_Ross_01022013/temp1/ICM23_A_index1_R1_001.fastq.gz
cat Sample_PRKI001A_index1/*_R2_*.fastq.gz > $rawData/Embryo.ICM/PE_100_fr.unstranded_Ross_01022013/temp1/ICM23_A_index1_R2_001.fastq.gz
cat Sample_PRKI001B_index3/*_R1_*.fastq.gz > $rawData/Embryo.ICM/PE_100_fr.unstranded_Ross_01022013/temp1/ICM24_B_index3_R1_001.fastq.gz
cat Sample_PRKI001B_index3/*_R2_*.fastq.gz > $rawData/Embryo.ICM/PE_100_fr.unstranded_Ross_01022013/temp1/ICM24_B_index3_R2_001.fastq.gz
cat Sample_PRKI001C_index4/*_R1_*.fastq.gz > $rawData/Embryo.ICM/PE_100_fr.unstranded_Ross_01022013/temp1/ICM31_C_index4_R1_001.fastq.gz
cat Sample_PRKI001C_index4/*_R2_*.fastq.gz > $rawData/Embryo.ICM/PE_100_fr.unstranded_Ross_01022013/temp1/ICM31_C_index4_R2_001.fastq.gz
bash ${script_path}/run_fastqc.sh "$rawData/Embryo.ICM/PE_100_fr.unstranded_Ross_01022013/temp1"
mkdir -p $prepData/Embryo.ICM/PE_100_fr.unstranded_Ross_01022013/fastq_data
mv $rawData/Embryo.ICM/PE_100_fr.unstranded_Ross_01022013/temp1/* $prepData/Embryo.ICM/PE_100_fr.unstranded_Ross_01022013/fastq_data/.

mkdir -p $rawData/Embryo.TE/PE_100_fr.unstranded_Ross_01022013/rawdata
cd $rawData/Embryo.TE/PE_100_fr.unstranded_Ross_01022013/rawdata
rsync -avzP -e ssh "tmansour@kaiden.genomecenter.ucdavis.edu:/finno/data/Ross_Data/Project_Ross_P/{Sample_PRKI001D_index5,Sample_PRKI001E_index6,Sample_PRKI001F_index7}" .
mkdir -p $rawData/Embryo.TE/PE_100_fr.unstranded_Ross_01022013/temp1
cat Sample_PRKI001D_index5/*_R1_*.fastq.gz > $rawData/Embryo.TE/PE_100_fr.unstranded_Ross_01022013/temp1/TE2313_D_index5_R1_001.fastq.gz
cat Sample_PRKI001D_index5/*_R2_*.fastq.gz > $rawData/Embryo.TE/PE_100_fr.unstranded_Ross_01022013/temp1/TE2313_D_index5_R2_001.fastq.gz
cat Sample_PRKI001E_index6/*_R1_*.fastq.gz > $rawData/Embryo.TE/PE_100_fr.unstranded_Ross_01022013/temp1/TE2413_E_index6_R1_001.fastq.gz
cat Sample_PRKI001E_index6/*_R2_*.fastq.gz > $rawData/Embryo.TE/PE_100_fr.unstranded_Ross_01022013/temp1/TE2413_E_index6_R2_001.fastq.gz
cat Sample_PRKI001F_index7/*_R1_*.fastq.gz > $rawData/Embryo.TE/PE_100_fr.unstranded_Ross_01022013/temp1/TE31_F_index7_R1_001.fastq.gz
cat Sample_PRKI001F_index7/*_R2_*.fastq.gz > $rawData/Embryo.TE/PE_100_fr.unstranded_Ross_01022013/temp1/TE31_F_index7_R2_001.fastq.gz
bash ${script_path}/run_fastqc.sh "$rawData/Embryo.TE/PE_100_fr.unstranded_Ross_01022013/temp1"
mkdir -p $prepData/Embryo.TE/PE_100_fr.unstranded_Ross_01022013/fastq_data
mv $rawData/Embryo.TE/PE_100_fr.unstranded_Ross_01022013/temp1/* $prepData/Embryo.TE/PE_100_fr.unstranded_Ross_01022013/fastq_data/.




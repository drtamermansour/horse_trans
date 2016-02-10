horse_trans="$1"

sudo apt-get update && \
sudo apt-get -y install screen git curl gcc make g++ python-dev unzip \
         default-jre pkg-config libncurses5-dev r-base-core r-cran-gplots \
         python-matplotlib python-pip python-virtualenv sysstat fastqc \
         trimmomatic bowtie samtools blast2

source $horse_trans/config.txt
echo -e ">PrefixPE/1\nTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n>PrefixPE/2\nGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n" >  $script_path/TruSeq3-PE.fa


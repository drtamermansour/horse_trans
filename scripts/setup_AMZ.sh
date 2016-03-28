horse_trans="$1"

sudo apt-get update && \
sudo apt-get -y install screen git curl gcc make g++ python-dev unzip \
         default-jre pkg-config libncurses5-dev r-base-core r-cran-gplots \
         python-matplotlib python-pip python-virtualenv sysstat fastqc \
         trimmomatic bowtie samtools blast2 

source $horse_trans/config.txt
echo -e ">PrefixPE/1\nTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n>PrefixPE/2\nGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n" >  $script_path/TruSeq3-PE.fa

sudo apt-get install bowtie2 bedtools tophat

sudo apt-get --purge remove bowtie2
cd /usr/src/
sudo wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.7/bowtie2-2.2.7-source.zip
sudo unzip bowtie2-2.2.7-source.zip
cd bowtie2-2.2.7
sudo make
sudo echo "export PATH=$PATH:/usr/src/bowtie2-2.2.7" >> ~/.bashrc
source ~/.bashrc

sudo apt-get --purge remove tophat
cd /usr/src/
sudo wget http://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.0.Linux_x86_64.tar.gz
sudo tar -zxvf tophat-2.1.0.Linux_x86_64.tar.gz
cd /usr/bin
sudo ln -s /usr/src/tophat-2.1.0.Linux_x86_64/tophat2 ./tophat



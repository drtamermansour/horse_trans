if (!require("ggplot2")) {
    install.packages("ggplot2", dependencies = TRUE, lib='~/R/v3.0.1/library', contriburl=contrib.url('http://cran.r-project.org/'))
    library(ggplot2, lib.loc = "~/R/v3.0.1/library")
}
if (!require("reshape2")) {
    install.packages("reshape2", dependencies = TRUE, lib='~/R/v3.0.1/library', contriburl=contrib.url('http://cran.r-project.org/'))
    library(reshape2, lib.loc = "~/R/v3.0.1/library")
}

args=(commandArgs(TRUE));
data_chr<-read.table(args[1], header=T)

##making the genes/chr/class code plot with chr in order
#Making the plot just against NCBI class codes
keeps_NCBI <- c("transcript.ID", "chr", "NCBI.class_code")
data_chr_NCBI <- data_chr[keeps_NCBI]
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28", "chr29", "chr30", "chr31", "chrUn", "chrX")
chrs_N <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "Un", "X")
#pdf("transcriptsCompareNCBIperChr.pdf")
#ggplot(subset(data_chr_NCBI,NCBI.class_code %in% c("j","=", "u"))) + scale_x_discrete(limits = chrs, labels = chrs_N) +
#  geom_bar(aes(chr, group=NCBI.class_code, colour=NCBI.class_code)) + theme(legend.position="top") + ylab("gene count") +
#  scale_colour_discrete(name  ="class code",
#                        labels=c("match","similar","novel","others"),
#                        expand=2) +
#  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
#  theme(legend.text = element_text(colour="black", size = 16)) +
#  theme(axis.text = element_text(colour="black", size = 12)) +
#  theme(axis.title = element_text(colour="black", size = 14))
#dev.off()

levels(data_chr_NCBI$NCBI.class_code) <- c(levels(data_chr_NCBI$NCBI.class_code), "others")
data_chr_NCBI$NCBI.class_code[!(data_chr_NCBI$NCBI.class_code == "j" | data_chr_NCBI$NCBI.class_code == "=" | data_chr_NCBI$NCBI.class_code == "u")]="others"
#pdf("transcriptsCompareNCBIperChr2.pdf")
#ggplot(data_chr_NCBI) + scale_x_discrete(limits = chrs, labels = chrs_N) +
#  geom_bar(aes(chr, group=NCBI.class_code, colour=NCBI.class_code)) + theme(legend.position="top") + ylab("gene count") +
#  scale_colour_discrete(name  ="class code",
#  labels=c("match","similar","novel","others"),
#  expand=2) +
#  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
#  theme(legend.text = element_text(colour="black", size = 16)) +
#  theme(axis.text = element_text(colour="black", size = 12)) +
#  theme(axis.title = element_text(colour="black", size = 14))
#dev.off()

# Making step line for chr size
chromSizes<-read.table(args[2], header=F, col.names=c("chr","size"))
#chromSizes<-read.table("equCab2.chrom.sizes", header=F, col.names=c("chr","size"))
pdf(args[3])
ggplot(subset(data_chr_NCBI,chr %in% chrs)) + scale_x_discrete(limits = chrs, labels = chrs_N) +
geom_bar(aes(chr, group=NCBI.class_code, colour=NCBI.class_code)) + theme(legend.position="top") +
ylab("gene count") +
scale_colour_discrete(name="class code", labels=c("match","similar","novel","others"), expand=2) +
theme(legend.title = element_text(colour="black", size=18, face="bold")) +
theme(legend.text = element_text(colour="black", size = 16)) +
theme(axis.text = element_text(colour="black", size = 10)) +
theme(axis.title = element_text(colour="black", size = 14)) +
geom_line(data=subset(chromSizes,chr %in% chrs), aes(x=chr, y=size/50000, group=1), colour="yellow")
#geom_point(data=chromSizes, aes(x=chr, y=size/10000), colour="blue", size = 3)
dev.off()

# Making a figure with class_code/chr plot and chr size trend line
#library(Rmisc)
#multiplot(c,d,cols=1)


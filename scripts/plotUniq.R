if (!require("ggplot2")) {
    install.packages("ggplot2", dependencies = TRUE, lib='~/R/v3.2.0/library', contriburl=contrib.url('http://cran.r-project.org/'))
    library(ggplot2, lib.loc = "~/R/v3.2.0/library")
}
if (!require("reshape2")) {
    install.packages("reshape2", dependencies = TRUE, lib='~/R/v3.0.1/library', contriburl=contrib.url('http://cran.r-project.org/'))
    library(reshape2, lib.loc = "~/R/v3.0.1/library")
}
if (!require("plyr")) {
    install.packages("plyr", dependencies = TRUE, lib='~/R/v3.0.1/library', contriburl=contrib.url('http://cran.r-project.org/'))
    library(plyr, lib.loc = "~/R/v3.0.1/library")
}
if (!require("dplyr")) {
    install.packages("dplyr", dependencies = TRUE, lib='~/R/v3.0.1/library', contriburl=contrib.url('http://cran.r-project.org/'))
    library(dplyr, lib.loc = "~/R/v3.0.1/library")
}

args=(commandArgs(TRUE));
data_0<-read.table(args[1])
data_5<-read.table(args[2])

data_changed_0 <- cbind(as.data.frame(data_0[1:7,]),as.data.frame(data_0[8:14,]),as.data.frame(data_0[15:21,]),as.data.frame(data_0[22:28,]),as.data.frame(data_0[29:35,]),as.data.frame(data_0[36:42,]),as.data.frame(data_0[43:49,]),as.data.frame(data_0[50:56,]))   
data_changed_0 <- sapply(data_changed_0, as.character)
colnames(data_changed_0) <- data_changed_0[1,]
data_changed_0 <- as.data.frame(data_changed_0[-1,])
rownames(data_changed_0) <- c("genes","isoforms","unique_gene","unique_isoforms","not_unique_genes","not_unique_isoforms")
data_0 <-as.data.frame(t(data_changed_0))
write.table(data_0, "tissue_geneVSisoforms_0.txt")

data_changed_5 <- cbind(as.data.frame(data_5[1:7,]),as.data.frame(data_5[8:14,]),as.data.frame(data_5[15:21,]),as.data.frame(data_5[22:28,]),as.data.frame(data_5[29:35,]),as.data.frame(data_5[36:42,]),as.data.frame(data_5[43:49,]),as.data.frame(data_5[50:56,]))   
data_changed_5 <- sapply(data_changed_5, as.character)
colnames(data_changed_5) <- data_changed_5[1,]
data_changed_5 <- as.data.frame(data_changed_5[-1,])
rownames(data_changed_5) <- c("genes","isoforms","unique_gene","unique_isoforms","not_unique_genes","not_unique_isoforms")
data_5 <-as.data.frame(t(data_changed_5))
write.table(data_5, "tissue_geneVSisoforms_5.txt")

data_0 <- read.table("tissue_geneVSisoforms_0.txt",stringsAsFactors=FALSE)
data_0$not_unique_genes <- data_0$not_unique_genes * -1
data_0$not_unique_isoforms <- data_0$not_unique_isoforms * -1

data_5 <- read.table("tissue_geneVSisoforms_5.txt",stringsAsFactors=FALSE)
data_5$not_unique_genes <- data_5$not_unique_genes * -1
data_5$not_unique_isoforms <- data_5$not_unique_isoforms * -1

#more isoform plots
Absent_Uisoforms_data_0 <- as.data.frame(data_0$not_unique_isoforms)
rownames(Absent_Uisoforms_data_0) <- rownames(data_0) 
#NUisoforms_data_0 <- as.data.frame(data_0$unique_isoforms)
#rownames(NUisoforms_data_0) <- rownames(data_0)

#Absent_Uisoforms_data_5 <- as.data.frame(data_5$not_unique_isoforms)
#rownames(Absent_Uisoforms_data_5) <- rownames(data_5) 
NUisoforms_data_5 <- as.data.frame(data_5$unique_isoforms)
rownames(NUisoforms_data_5) <- rownames(data_5)

#This is figure 3B without isform count line
pdf(args[3],paper="a4r",width=10, height=5)
ggplot() +
  geom_bar(data=NUisoforms_data_5, aes(x=rownames(data_5),y=data_5$unique_isoforms,color="aliceblue"), stat="identity", position = "identity") +
  geom_bar(data=Absent_Uisoforms_data_0, aes(x=rownames(data_0),y=data_0$not_unique_isoforms,color="red"), stat="identity", position = "identity") + 
  ylab("Number of isoforms") + 
  scale_color_discrete(name="Unique isoforms", labels=c("present","absent")) +
  xlab("Tissue") +
  theme(legend.title = element_text(colour="black", size=14, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(colour="black", size = 9)) +
  theme(axis.title = element_text(colour="black", size = 14))
dev.off()

#Figure with isoform count trend line below:
currentPath=getwd()
setwd(args[4])
#load in uniquely present genes in each tissue
BrainStem_p <- read.table("BrainStem.isoform.expressed_uniqely_cutoff.5", stringsAsFactors=F)
Cerebellum_p <- read.table("Cerebellum.isoform.expressed_uniqely_cutoff.5", stringsAsFactors=F)
Muscle_p <- read.table("Muscle.isoform.expressed_uniqely_cutoff.5", stringsAsFactors=F)
Retina_p <- read.table("Retina.isoform.expressed_uniqely_cutoff.5", stringsAsFactors=F)
Embryo.ICM_p <- read.table("Embryo.ICM.isoform.expressed_uniqely_cutoff.5", stringsAsFactors=F)
Embryo.TE_p <- read.table("Embryo.TE.isoform.expressed_uniqely_cutoff.5", stringsAsFactors=F)
Skin_p <- read.table("Skin.isoform.expressed_uniqely_cutoff.5", stringsAsFactors=F)
Spinalcord_p <- read.table("SpinalCord.isoform.expressed_uniqely_cutoff.5", stringsAsFactors=F)

#load uniquely absent genes for each tissue
BrainStem_a <- read.table("BrainStem.isoform.notExpressed_uniqely_cutoff.0", stringsAsFactors=F)
Cerebellum_a <- read.table("Cerebellum.isoform.notExpressed_uniqely_cutoff.0", stringsAsFactors=F)
Muscle_a <- read.table("Muscle.isoform.notExpressed_uniqely_cutoff.0", stringsAsFactors=F)
Retina_a <- read.table("Retina.isoform.notExpressed_uniqely_cutoff.0", stringsAsFactors=F)
Embyro.ICM_a <- read.table("Embryo.ICM.isoform.notExpressed_uniqely_cutoff.0", stringsAsFactors=F)
Embryo.TE_a <- read.table("Embryo.TE.isoform.notExpressed_uniqely_cutoff.0", stringsAsFactors=F)
Skin_a <- read.table("Skin.isoform.notExpressed_uniqely_cutoff.0", stringsAsFactors=F)
Spinalcord_a <- read.table("SpinalCord.isoform.notExpressed_uniqely_cutoff.0", stringsAsFactors=F)

###Concatanate all expression data and add tissue label for downstream tracking
#label all with tissue
BrainStem_p_target=cbind(BrainStem_p[,c(1,2)],"BrainStem")
names(BrainStem_p_target)=c("transcript_name","TPM","Tissue")
Cerebellum_p_target=cbind(Cerebellum_p[,c(1,3)],"Cerebellum")
names(Cerebellum_p_target)=c("transcript_name","TPM","Tissue")
Embryo.ICM_p_target=cbind(Embryo.ICM_p[,c(1,4)],"Embryo.ICM")
names(Embryo.ICM_p_target)=c("transcript_name","TPM","Tissue")
Embryo.TE_p_target=cbind(Embryo.TE_p[,c(1,5)],"Embryo.TE")
names(Embryo.TE_p_target)=c("transcript_name","TPM","Tissue")
Muscle_p_target=cbind(Muscle_p[,c(1,6)],"Muscle")
names(Muscle_p_target)=c("transcript_name","TPM","Tissue")
Retina_p_target=cbind(Retina_p[,c(1,7)],"Retina")
names(Retina_p_target)=c("transcript_name","TPM","Tissue")
Skin_p_target=cbind(Skin_p[,c(1,8)],"Skin")
names(Skin_p_target)=c("transcript_name","TPM","Tissue")
Spinalcord_p_target=cbind(Spinalcord_p[,c(1,6)],"SpinalCord")
names(Spinalcord_p_target)=c("transcript_name","TPM","Tissue")
#rbind with label minus 1st row of column headers
all_unique <- rbind(BrainStem_p_target,Cerebellum_p_target,Embryo.ICM_p_target,Embryo.TE_p_target,
                    Muscle_p_target,Retina_p_target,Skin_p_target,Spinalcord_p_target)
#Melt data for manipulation
#melted <- melt(all_unique, id.vars=c("isoforms","gene.ID","exons","NCBI.ref_gene_id","chr","label","chr_start","chr_end"))
#Calculate the sum(TPM) and STDEV of each gene per tissue
#melted_tissue <- ddply(melted, c("variable"), summarise, sum = sum(value))
melted_tissue <- ddply(all_unique, c("Tissue"), summarise, sum = sum(TPM))

##unique present and solely absent isoforms with TPM line...FIGURE 3B
setwd(currentPath)
pdf(paste("withCounts",args[3],sep="."),paper="a4r",width=10, height=5)
ggplot() +
  geom_bar(data=NUisoforms_data_5, aes(x=rownames(data_5),y=data_5$unique_isoforms,color="aliceblue"), stat="identity", position = "identity") +
  geom_bar(data=Absent_Uisoforms_data_0, aes(x=rownames(data_0),y=data_0$not_unique_isoforms,color="red"), stat="identity", position = "identity") + 
  ylab("Number of isoforms") +
  scale_color_discrete(name="Unique isoforms", labels=c("present","absent"),expand = c(0,0)) +
  xlab("Tissue") +
  theme(legend.title = element_text(colour="black", size=14, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(colour="black", size = 9)) +
  theme(axis.title = element_text(colour="black", size = 14)) +
  geom_line(data=melted_tissue, aes(x=Tissue,y=sum / 50, group=1),colour="green") +
  scale_x_discrete(limits=c("BrainStem","Cerebellum","Embryo.ICM", "Embryo.TE","Muscle","Retina","Skin","SpinalCord"))
dev.off()



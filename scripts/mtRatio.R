if (!require("ggplot2")) {
    install.packages("ggplot2", dependencies = TRUE, lib='~/R/v3.2.0/library', contriburl=contrib.url('http://cran.r-project.org/'))
    library(ggplot2, lib.loc = "~/R/v3.2.0/library")
}
if (!require("reshape2")) {
    install.packages("reshape2", dependencies = TRUE, lib='~/R/v3.2.0/library', contriburl=contrib.url('http://cran.r-project.org/'))
    library(reshape2, lib.loc = "~/R/v3.2.0/library")
}
if (!require("dplyr")) {
    install.packages("dplyr", dependencies = TRUE, lib='~/R/v3.2.0/library', contriburl=contrib.url('http://cran.r-project.org/'))
    library(dplyr, lib.loc = "~/R/v3.2.0/library")
}

#read in dataset with expression values for each term and a dataset with annotation for each gene
args=(commandArgs(TRUE));
tissues <- read.table(args[1],header=T)
colnames(tissues)=gsub("Embryo.", "Emb.", colnames(tissues))

mt <- tissues[1:2,2:9]
mt_TPM <- as.data.frame(colSums(mt))
nuclear <- tissues [-(1:2),2:9]
nuclear_TPM <- as.data.frame(colSums(nuclear))

t_out <- as.data.frame(cbind(mt_TPM,nuclear_TPM))
colnames(t_out) <- c("mitochondrial", "nuclear")

#reshaping data for geom_bar
df.melt <- melt(as.matrix(t_out))
bar <- group_by(df.melt, Var2, Var1)

#plotting nuclear vs mitochondrial origin transcripts per tissue as proportion
png(filename=args[2], width=900, height=600)
ggplot(bar, aes(x=Var1, y=value / 1000000, fill=factor(Var2)))+
  geom_bar(position="stack", stat="identity") + 
  ylab("Proportion of transcriptional output") + xlab("Tissue") +
  guides(fill=guide_legend(title="gene origin",
                           labels=c("nuclear", "mitochondrial"))) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(legend.title = element_text(colour="black", size = 18)) +
  theme(axis.text = element_text(colour="black", size = 16)) +
  theme(axis.title = element_text(colour="black", size = 20))
graphics.off()  # close the PNG device


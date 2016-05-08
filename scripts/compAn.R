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
#data_chr<-read.table("nonGuided_Cufflinks.nonGuided_Cuffmerge.merge.reduced", header=T)
keeps <- c("transcript.ID", "NCBI.class_code", "ensGTF_file.class_code", "ISME.PBMC.class_code", "Hestand_2014.class_code")
selected <- data_chr[keeps]
selected_class <- selected[ ,2:5]
colnames(selected_class) <- c("NCBI","ENSEMBL","ISME","Hestand")
flip<-t(selected_class)
df.melt <- melt(flip, id=rownames(flip))
levels(df.melt$value) <- c(levels(df.melt$value), "others")
df.melt$value[!(df.melt$value == "j" | df.melt$value == "=" | df.melt$value == "u")]="others"
#The bar graph comparing dbs (Venn substitution)
pdf(args[2])
ggplot(df.melt,aes(Var1, fill=value)) + 
  geom_bar(position="dodge") + xlab("database") + ylab("transcript count") +
  scale_fill_manual(values=c("tomato","green3","dodgerblue1","gray"),
                    name="class code",
                    breaks=c("=", "j", "u", "others"),
                    labels=c("match", "similar","novel", "others"))
dev.off()


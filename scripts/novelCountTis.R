if (!require("RColorBrewer")) {
    install.packages("RColorBrewer", dependencies = TRUE, lib='~/R/v3.2.0/library', contriburl=contrib.url('http://cran.r-project.org/'))
    library(RColorBrewer, lib.loc = "~/R/v3.2.0/library")
}
if (!require("ggplot2")) {
    install.packages("ggplot2", dependencies = TRUE, lib='~/R/v3.2.0/library', contriburl=contrib.url('http://cran.r-project.org/'))
    library(ggplot2, lib.loc = "~/R/v3.2.0/library")
}
if (!require("reshape2")) {
    install.packages("reshape2", dependencies = TRUE, lib='~/R/v3.2.0/library', contriburl=contrib.url('http://cran.r-project.org/'))
    library(reshape2, lib.loc = "~/R/v3.2.0/library")
}
if (!require("plyr")) {
    install.packages("plyr", dependencies = TRUE, lib='~/R/v3.2.0/library', contriburl=contrib.url('http://cran.r-project.org/'))
    library(plyr, lib.loc = "~/R/v3.2.0/library")
}

args=(commandArgs(TRUE));
reference_sup <- read.table(args[1],stringsAsFactors=FALSE, header=T)
reference_cons <- read.table(args[2],stringsAsFactors=FALSE, header=T)
reference_unsup <- read.table(args[3],stringsAsFactors=FALSE, header=T)

keeps_exons <- c("exons","transcript.ID")
novel_sup_exons <- reference_sup[keeps_exons]
novel_cons_exons <- reference_cons[keeps_exons]
novel_unsup_exons <- reference_unsup[keeps_exons]

### Extracting expression information about the novel isoforms
#start by uploading TPM values for all isoforms
tissues <- read.table(args[4], header=T)
#merge the TPM values with the novel gene by "geneName"
m_express_sup <- merge(novel_sup_exons,tissues,by.x="transcript.ID",by.y="isoformName")
m_express_cons <- merge(novel_cons_exons,tissues,by.x="transcript.ID",by.y="isoformName")
m_express_unsup <- merge(novel_unsup_exons,tissues,by.x="transcript.ID",by.y="isoformName")

#reshaping data for geom_bar
df.melt_sup <- melt(m_express_sup, id=c("transcript.ID","exons"))
df.melt_cons <- melt(m_express_cons, id=c("transcript.ID","exons"))
df.melt_unsup <- melt(m_express_unsup, id=c("transcript.ID","exons"))

#combine and label by novel gene category
df.melt_sup$cat <- 'Category I'
df.melt_cons$cat <- 'Category II'
df.melt_unsup$cat <- 'Category III'
df.melt <- rbind(df.melt_sup,df.melt_cons,df.melt_unsup)

## Plot sum of TPM of novel genes with exon number subset per tissue
my_pal <- brewer.pal(6,"Set3")
#getting gene numbers
geneN <- subset(df.melt, value > 0.5)
#geneN <- geneN[!duplicated(geneN),]
count_gene <- count(geneN,c("variable","cat","exons"))
# update a new version of merged data by changing all exon counts to be 6 for easier visualiztion 
df.melt2=df.melt
df.melt2$exons[which(df.melt2$exons > 5)]=6 
# plot
png(filename=args[5], width=1000, height=750)
ggplot(df.melt2, aes(variable,group=exons,fill=exons)) +
  geom_bar(aes(weight=value),position="stack") + xlab("Tissue") +
  ylab("cumulative expression (TPM)") +
 # scale_fill_gradientn(breaks=c(1,2,3,4,5,6), labels=c("1","2","3","4","5",">5"), colours=my_pal, guide="legend") +
  scale_fill_gradientn(breaks=c(6,5,4,3,2,1), labels=c(">5","5","4","3","2","1"), colours=my_pal, guide="legend") +
  guides(fill=guide_legend(title="exon number")) +
  #geom_point(data=count_gene, aes(x=variable, y=freq * 50,size=exons, group=1)) +
  #scale_size_continuous(limits=c(1,2,3,4,5,6),
                      #labels=c("1","2","3", "4","5",">5")) +
  facet_grid(~ cat) +
  theme(axis.text.x = element_text(colour="black", size =10,angle=90 )) +
  theme(axis.title = element_text(colour="black", size = 14))
graphics.off()  # close the PNG device

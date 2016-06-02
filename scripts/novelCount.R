if (!require("RColorBrewer")) {
    install.packages("RColorBrewer", dependencies = TRUE, lib='~/R/v3.2.0/library', contriburl=contrib.url('http://cran.r-project.org/'))
    library(RColorBrewer, lib.loc = "~/R/v3.2.0/library")
}
if (!require("ggplot2")) {
    install.packages("ggplot2", dependencies = TRUE, lib='~/R/v3.2.0/library', contriburl=contrib.url('http://cran.r-project.org/'))
    library(ggplot2, lib.loc = "~/R/v3.2.0/library")
}
args=(commandArgs(TRUE));
reference_sup <- read.table(args[1],stringsAsFactors=FALSE, header=T)
reference_cons <- read.table(args[2],stringsAsFactors=FALSE, header=T)
reference_unsup <- read.table(args[3],stringsAsFactors=FALSE, header=T)

keeps_exons <- c("exons","transcript.ID")
novel_sup_exons <- reference_sup[keeps_exons]
novel_sup_exons["label"]<- "novel_supported"
novel_cons_exons <- reference_cons[keeps_exons]
novel_cons_exons["label"]<- "novel_conserved"
novel_unsup_exons <- reference_unsup[keeps_exons]
novel_unsup_exons["label"]<- "novel_unsupported"

all_novel_exons=rbind(novel_sup_exons,novel_cons_exons,novel_unsup_exons)
#making the figure with ability to see all exons
colourCount = length(unique(all_novel_exons$exons))
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
png(filename=args[4], width=800, height=750)
ggplot(all_novel_exons) + 
  stat_count(aes(label,fill=factor(exons))) + ylab("Transcript count") + xlab("Novel gene category") +
  scale_y_discrete("Transcript count", limits=(seq(0,15000,1000))) +
  guides(fill=guide_legend(title="Number of exons")) +
  theme(axis.text.y = element_text(colour="black", size = 14, angle = 90)) +
  theme(axis.text.x = element_text(colour="black", size = 12)) +
  theme(axis.text.y = element_text(colour="black", size = 10, angle=0)) +
  theme(axis.title = element_text(colour="black", size = 14)) +
  theme(legend.text = element_text(colour="black", size = 10)) +
  scale_x_discrete("Novel gene category", labels = c("novel_conserved" = "Category II","novel_supported" = " Category I",
                                                     "novel_unsupported" = "Category III"),
                   limits=c("novel_supported","novel_conserved","novel_unsupported")) +
  scale_fill_manual(breaks=c(1,2,3,4,5,10,28),
                    labels=c("1","2","3","4","5","6-10","11-28"),
                    values = getPalette(colourCount))
graphics.off()  # close the PNG device

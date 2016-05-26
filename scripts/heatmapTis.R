##Figure 3A: Tissue specificity, heatmap of genes
###Tissue specific heatmap
if (!require("reshape2")) {
    install.packages("reshape2", dependencies = TRUE, lib='~/R/v3.2.0/library', contriburl=contrib.url('http://cran.r-project.org/'))
    library(reshape2, lib.loc = "~/R/v3.2.0/library")
}
if (!require("plyr")) {
    install.packages("plyr", dependencies = TRUE, lib='~/R/v3.2.0/library', contriburl=contrib.url('http://cran.r-project.org/'))
    library(plyr, lib.loc = "~/R/v3.2.0/library")
}

#Read in data
args=(commandArgs(TRUE));
data <- read.table(args[1], header=T)
#Melt data for manipulation
melted <- melt(data, id.vars=c("geneName"))
#Calculate the sum(TPM) and STDEV of each gene per tissue
melted_new<- ddply(melted, c("geneName"), summarise,
                   sum = sum(value), sd = sd(value))
#Add the column of sum and sd to the original TPM values table
complete <- merge(data,melted_new,by="geneName")
#Subset data based on if sum>50 and sd>50
sub <- subset(complete, c(sum > 200 & sd > 200))
#make row.names the geneName
rownames(sub)<-sub$geneName
# making the matrix for the heatmap
#disable scientific notation so no "e+/-"
options("scipen"=100, "digits"=4)
datanumbers <- data.matrix(sub[,2:9])

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("Blue", "white", "Red"))(n = 18)
#making the heatmap
if (!require("gplots")) {
    install.packages("gplots", dependencies = TRUE, lib='~/R/v3.2.0/library', contriburl=contrib.url('http://cran.r-project.org/'))
    library(gplots, lib.loc = "~/R/v3.2.0/library")
}
if (!require("RColorBrewer")) {
    install.packages("RColorBrewer", dependencies = TRUE, lib='~/R/v3.2.0/library', contriburl=contrib.url('http://cran.r-project.org/'))
    library(RColorBrewer, lib.loc = "~/R/v3.2.0/library")
}
#if (!require("svDialogs")) {
#    install.packages("svDialogs", dependencies = TRUE, lib='~/R/v3.2.0/library', contriburl=contrib.url('http://cran.r-project.org/'))
#    library(svDialogs, lib.loc = "~/R/v3.2.0/library")
#}

#Manually cluster to your liking...allows use you to choose which correlations work best 
#with your data depending on linear or monotonic relationship of variables
# Row clustering(genes)...pearson seems to work best here 
hr <- hclust(as.dist(1-cor(t(datanumbers), method="pearson")), method="average")
# Column clustering(tissues)...spearman seems to work best here 
hc <- hclust(as.dist(1-cor(datanumbers, method="spearman")), method="average")

## Plot heatmap
par(mar=c(7,4,4,2)+0.1) 
png(filename=args[2], width=800, height=750)
col_breaks <- c(1:10,20,30,40,50,60,70,80,90,100)
heatmap.2(datanumbers,    # data matrix
          #cellnote = mat_data,  # same data set for cell labels
          #main = "Rank", # heat map title
          #notecex=0.4,
          #cexRow=0.6,
          cexCol=2,
          scale="none",
          #notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,12),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="column",     # only draw a row dendrogram
          Colv=as.dendrogram(hc),
          Rowv=as.dendrogram(hr),
          hclustfun = hclust,
          labRow = NULL,
          key.xlab = "TPM",
          key.title = NULL)            # turn off column clustering
graphics.off()  # close the PNG device
## Return matrix with row/column sorting as in heatmap
write.csv(datanumbers[rev(hr$labels[hr$order]), hc$labels[hc$order]],args[3])
#match the XLOC_* names with real gene names and add other annotation information
#So you can see gene name and chr location,strand, and how it compares in other databases
hmap_order <- read.csv(args[3], header=T)
annotated <- read.table(args[4], header=T)
m <- merge(hmap_order, annotated, by.x="X",by.y="gene.ID", sort=F)
write.csv(m,paste("ann",args[3],sep="."))

